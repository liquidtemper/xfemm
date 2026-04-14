from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
from typing import Dict, Iterable, List, Tuple

from generate_virtual_gap_cases import main as generate_cases

try:
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"matplotlib is required for benchmark plots: {exc}")


ROOT = Path(__file__).resolve().parent
CASE_DIR = ROOT / "cases"
OUT_DIR = ROOT / "out"
FMESHER = Path("C:/GIT/xfemm/cfemm/bin/fmesher")
FSOLVER = Path("C:/GIT/xfemm/cfemm/bin/fsolver")
FPPROC_VGAP = Path("C:/GIT/xfemm/cfemm/bin/fpproc-vgap")


CASE_SPECS = {
    "case1": {
        "variants": [
            "case1_steel_steel_nominal.fem",
            "case1_steel_steel_vgap_100um.fem",
            "case1_steel_steel_truegap_100um.fem",
        ],
        "field_box": (-5.0, -5.0, 85.0, 105.0),
        "field_grid": (120, 120),
        "probe_points": [(70.0, 50.0)],
        "metrics": {
            "circuit": "coil",
            "line_bn_avg": ((61.0, 50.0), (79.0, 50.0)),
        },
    },
    "case2": {
        "variants": [
            "case2_surface_layer_nominal.fem",
            "case2_surface_layer_vgap_100um.fem",
        ],
        "field_box": (10.0, 0.0, 90.0, 95.0),
        "field_grid": (120, 120),
        "probe_points": [(50.0, 70.1), (50.0, 70.6)],
        "metrics": {
            "circuit": "coil",
            "line_bn_avg": ((40.0, 70.1), (60.0, 70.1)),
            "peak_window": (38.0, 70.2, 62.0, 71.0),
        },
    },
    "case3": {
        "variants": [
            "case3_pm_interface_nominal.fem",
            "case3_pm_interface_vgap_100um.fem",
        ],
        "field_box": (15.0, -5.0, 65.0, 55.0),
        "field_grid": (120, 120),
        "probe_points": [(40.0, 44.0), (52.0, 25.0)],
        "metrics": {
            "line_bn_avg": ((40.0, 38.0), (40.0, 50.0)),
        },
    },
}


def win_to_wsl(path: Path) -> str:
    path = path.resolve()
    drive = path.drive.rstrip(":").lower()
    tail = path.as_posix().split(":/", 1)[1]
    return f"/mnt/{drive}/{tail}"


def run_case(fem_file: Path, spec: Dict, work_dir: Path) -> Tuple[Dict[str, float], List[Tuple[float, float, float]]]:
    metrics_csv = work_dir / f"{fem_file.stem}_metrics.csv"
    field_csv = work_dir / f"{fem_file.stem}_field.csv"
    fem_out = work_dir / f"{fem_file.stem}.fem"
    subprocess.run(
        ["wsl.exe", "bash", "-lc", f'cp "{win_to_wsl(fem_file)}" "{win_to_wsl(fem_out)}"'],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    fem_stem = win_to_wsl(fem_out.with_suffix(""))
    subprocess.run(
        ["wsl.exe", "bash", "-lc", f'"{win_to_wsl(FMESHER)}" "{win_to_wsl(fem_out)}"'],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    subprocess.run(
        ["wsl.exe", "bash", "-lc", f'"{win_to_wsl(FSOLVER)}" "{fem_stem}"'],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    case_name = fem_file.stem.split("_", 1)[0]
    ans_path = fem_out.with_suffix(".ans")
    subprocess.run(
        [
            "wsl.exe",
            "bash",
            "-lc",
            f'"{win_to_wsl(FPPROC_VGAP)}" "{case_name}" "{win_to_wsl(ans_path)}" '
            f'"{win_to_wsl(metrics_csv)}" "{win_to_wsl(field_csv)}"',
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    metrics: Dict[str, float] = {}
    with metrics_csv.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            metrics[row["metric"]] = float(row["value"])

    field_rows: List[Tuple[float, float, float]] = []
    with field_csv.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            field_rows.append((float(row["x"]), float(row["y"]), float(row["bmag"])))
    return metrics, field_rows


def save_plot(case_name: str, variants: List[str], fields: Dict[str, List[Tuple[float, float, float]]], spec: Dict) -> None:
    nx, ny = spec["field_grid"]
    fig, axes = plt.subplots(1, len(variants), figsize=(5 * len(variants), 4), constrained_layout=True)
    if len(variants) == 1:
        axes = [axes]
    vmin = min(val for rows in fields.values() for _, _, val in rows)
    vmax = max(val for rows in fields.values() for _, _, val in rows)
    xmin, ymin, xmax, ymax = spec["field_box"]
    for ax, variant in zip(axes, variants):
        rows = fields[variant]
        grid = [[0.0] * nx for _ in range(ny)]
        for idx, (_, _, bmag) in enumerate(rows):
            iy = idx // nx
            ix = idx % nx
            grid[iy][ix] = bmag
        im = ax.imshow(
            grid,
            origin="lower",
            extent=(xmin, xmax, ymin, ymax),
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
            cmap="viridis",
        )
        ax.set_title(variant.replace(".fem", ""))
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
    fig.colorbar(im, ax=axes, label="|B| [T]")
    fig.savefig(OUT_DIR / f"{case_name}_field_compare.png", dpi=180)
    plt.close(fig)


def main() -> None:
    generate_cases()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows: List[Dict[str, object]] = []
    all_metrics: Dict[str, Dict[str, float]] = {}

    for case_name, spec in CASE_SPECS.items():
        fields: Dict[str, List[Tuple[float, float, float]]] = {}
        for variant in spec["variants"]:
            fem_file = CASE_DIR / variant
            metrics, field_rows = run_case(fem_file, spec, OUT_DIR)
            fields[variant] = field_rows
            all_metrics[variant] = metrics
            row: Dict[str, object] = {"case": case_name, "variant": variant}
            row.update(metrics)
            summary_rows.append(row)
        save_plot(case_name, spec["variants"], fields, spec)

    case1_vgap = all_metrics["case1_steel_steel_vgap_100um.fem"]
    case1_true = all_metrics["case1_steel_steel_truegap_100um.fem"]
    if "inductance" in case1_vgap and "inductance" in case1_true and case1_true["inductance"] != 0:
        err = 100.0 * abs(case1_vgap["inductance"] - case1_true["inductance"]) / case1_true["inductance"]
        summary_rows.append({
            "case": "case1",
            "variant": "vgap_vs_truegap_error",
            "inductance_percent_error": err,
        })

    summary_csv = OUT_DIR / "virtual_gap_benchmark_summary.csv"
    keys = sorted({k for row in summary_rows for k in row.keys()})
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(summary_rows)

    with (OUT_DIR / "virtual_gap_benchmark_summary.json").open("w") as f:
        json.dump(summary_rows, f, indent=2)


if __name__ == "__main__":
    main()
