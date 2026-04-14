from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import math


@dataclass
class BoundaryProp:
    name: str
    bdry_type: int = 0
    a0: float = 0.0
    a1: float = 0.0
    a2: float = 0.0
    phi: float = 0.0
    c0: float = 0.0
    c0i: float = 0.0
    c1: float = 0.0
    c1i: float = 0.0
    mu_ssd: float = 0.0
    sigma_ssd: float = 0.0


@dataclass
class Material:
    name: str
    mu_x: float
    mu_y: float
    h_c: float = 0.0
    h_c_angle: float = 0.0
    j_re: float = 0.0
    j_im: float = 0.0
    sigma: float = 0.0
    d_lam: float = 0.0
    phi_h: float = 0.0
    phi_hx: float = 0.0
    phi_hy: float = 0.0
    lam_type: int = 0
    lam_fill: float = 1.0
    n_strands: int = 0
    wire_d: float = 0.0
    bh_points: Sequence[Tuple[float, float]] = field(default_factory=list)


@dataclass
class Circuit:
    name: str
    amps_re: float
    amps_im: float = 0.0
    circuit_type: int = 1


@dataclass
class BlockLabel:
    x: float
    y: float
    material: Optional[str]
    meshsize: float = -1.0
    circuit: Optional[str] = None
    magdir: float = 0.0
    group: int = 0
    turns: int = 1
    is_external: bool = False
    is_default: bool = False


@dataclass
class Segment:
    n0: int
    n1: int
    max_side_length: float = -1.0
    boundary: Optional[str] = None
    hidden: bool = False
    group: int = 0
    virtual_gap_enabled: bool = False
    virtual_gap_physical: float = 0.0
    virtual_gap_mesh: float = -1.0


class FemModel:
    def __init__(self, *, comment: str, depth: float = 20.0, units: str = "millimeters",
                 precision: float = 1e-8, min_angle: float = 30.0) -> None:
        self.comment = comment
        self.depth = depth
        self.units = units
        self.precision = precision
        self.min_angle = min_angle
        self.boundaries: List[BoundaryProp] = []
        self.materials: List[Material] = []
        self.circuits: List[Circuit] = []
        self.labels: List[BlockLabel] = []
        self.segments: List[Segment] = []
        self._point_index: Dict[Tuple[int, int], int] = {}
        self.points: List[Tuple[float, float]] = []

    def add_point(self, x: float, y: float) -> int:
        key = (round(x * 1_000_000), round(y * 1_000_000))
        idx = self._point_index.get(key)
        if idx is not None:
            return idx
        idx = len(self.points)
        self._point_index[key] = idx
        self.points.append((x, y))
        return idx

    def add_segment(self, p0: Tuple[float, float], p1: Tuple[float, float], **kwargs) -> Segment:
        n0 = self.add_point(*p0)
        n1 = self.add_point(*p1)
        seg = Segment(n0=n0, n1=n1, **kwargs)
        self.segments.append(seg)
        return seg

    def add_polyline(self, pts: Sequence[Tuple[float, float]], close: bool = True, **kwargs) -> None:
        for i in range(len(pts) - 1):
            self.add_segment(pts[i], pts[i + 1], **kwargs)
        if close and len(pts) > 2:
            self.add_segment(pts[-1], pts[0], **kwargs)

    def add_material(self, material: Material) -> None:
        self.materials.append(material)

    def add_boundary(self, boundary: BoundaryProp) -> None:
        self.boundaries.append(boundary)

    def add_circuit(self, circuit: Circuit) -> None:
        self.circuits.append(circuit)

    def add_label(self, label: BlockLabel) -> None:
        self.labels.append(label)

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        boundary_map = {b.name: i + 1 for i, b in enumerate(self.boundaries)}
        material_map = {m.name: i + 1 for i, m in enumerate(self.materials)}
        circuit_map = {c.name: i + 1 for i, c in enumerate(self.circuits)}

        def fmt(value: float) -> str:
            return f"{value:.17g}"

        def ext_default(lbl: BlockLabel) -> int:
            out = 0
            if lbl.is_external:
                out |= 0x01
            if lbl.is_default:
                out |= 0x02
            return out

        with path.open("w", newline="\n") as f:
            f.write("[Format]      =  4.0\n")
            f.write("[Frequency]   =  0\n")
            f.write(f"[Precision]   =  {fmt(self.precision)}\n")
            f.write(f"[MinAngle]    =  {fmt(self.min_angle)}\n")
            f.write(f"[Depth]       =  {fmt(self.depth)}\n")
            f.write(f"[LengthUnits] =  {self.units}\n")
            f.write("[ProblemType] =  planar\n")
            f.write("[Coordinates] =  cartesian\n")
            f.write("[ACSolver]    =  0\n")
            f.write('[PrevSoln]    = ""\n')
            f.write("[PrevType]    =  0\n")
            f.write(f'[Comment]     =  "{self.comment}"\n')
            f.write("[PointProps]  =  0\n")
            f.write(f"[BdryProps]   = {len(self.boundaries)}\n")
            for b in self.boundaries:
                f.write("  <BeginBdry>\n")
                f.write(f'    <BdryName> = "{b.name}"\n')
                f.write(f"    <BdryType> = {b.bdry_type}\n")
                f.write(f"    <A_0> = {fmt(b.a0)}\n")
                f.write(f"    <A_1> = {fmt(b.a1)}\n")
                f.write(f"    <A_2> = {fmt(b.a2)}\n")
                f.write(f"    <Phi> = {fmt(b.phi)}\n")
                f.write(f"    <c0> = {fmt(b.c0)}\n")
                f.write(f"    <c0i> = {fmt(b.c0i)}\n")
                f.write(f"    <c1> = {fmt(b.c1)}\n")
                f.write(f"    <c1i> = {fmt(b.c1i)}\n")
                f.write(f"    <Mu_ssd> = {fmt(b.mu_ssd)}\n")
                f.write(f"    <Sigma_ssd> = {fmt(b.sigma_ssd)}\n")
                f.write("  <EndBdry>\n")
            f.write(f"[BlockProps]  = {len(self.materials)}\n")
            for m in self.materials:
                f.write("  <BeginBlock>\n")
                f.write(f'    <BlockName> = "{m.name}"\n')
                f.write(f"    <Mu_x> = {fmt(m.mu_x)}\n")
                f.write(f"    <Mu_y> = {fmt(m.mu_y)}\n")
                f.write(f"    <H_c> = {fmt(m.h_c)}\n")
                f.write(f"    <H_cAngle> = {fmt(m.h_c_angle)}\n")
                f.write(f"    <J_re> = {fmt(m.j_re)}\n")
                f.write(f"    <J_im> = {fmt(m.j_im)}\n")
                f.write(f"    <Sigma> = {fmt(m.sigma)}\n")
                f.write(f"    <d_lam> = {fmt(m.d_lam)}\n")
                f.write(f"    <Phi_h> = {fmt(m.phi_h)}\n")
                f.write(f"    <Phi_hx> = {fmt(m.phi_hx)}\n")
                f.write(f"    <Phi_hy> = {fmt(m.phi_hy)}\n")
                f.write(f"    <LamType> = {m.lam_type}\n")
                f.write(f"    <LamFill> = {fmt(m.lam_fill)}\n")
                f.write(f"    <NStrands> = {m.n_strands}\n")
                f.write(f"    <WireD> = {fmt(m.wire_d)}\n")
                f.write(f"    <BHPoints> = {len(m.bh_points)}\n")
                for b_val, h_val in m.bh_points:
                    f.write(f"      {fmt(b_val)}\t{fmt(h_val)}\n")
                f.write("  <EndBlock>\n")
            f.write(f"[CircuitProps]  = {len(self.circuits)}\n")
            for c in self.circuits:
                f.write("  <BeginCircuit>\n")
                f.write(f'    <CircuitName> = "{c.name}"\n')
                f.write(f"    <TotalAmps_re> = {fmt(c.amps_re)}\n")
                f.write(f"    <TotalAmps_im> = {fmt(c.amps_im)}\n")
                f.write(f"    <CircuitType> = {c.circuit_type}\n")
                f.write("  <EndCircuit>\n")
            f.write(f"[NumPoints] = {len(self.points)}\n")
            for x, y in self.points:
                f.write(f"{fmt(x)}\t{fmt(y)}\t0\t0\n")
            f.write(f"[NumSegments] = {len(self.segments)}\n")
            for seg in self.segments:
                bdry = boundary_map.get(seg.boundary, 0)
                f.write(
                    f"{seg.n0}\t{seg.n1}\t{fmt(seg.max_side_length)}\t{bdry}\t"
                    f"{1 if seg.hidden else 0}\t{seg.group}\t"
                    f"{1 if seg.virtual_gap_enabled else 0}\t"
                    f"{fmt(seg.virtual_gap_physical)}\t{fmt(seg.virtual_gap_mesh)}\n"
                )
            f.write("[NumArcSegments] = 0\n")
            f.write("[NumHoles] = 0\n")
            f.write(f"[NumBlockLabels] = {len(self.labels)}\n")
            for lbl in self.labels:
                block = 0 if lbl.material is None else material_map[lbl.material]
                meshsize = -1.0 if lbl.meshsize <= 0 else lbl.meshsize
                incircuit = 0 if lbl.circuit is None else circuit_map[lbl.circuit]
                f.write(
                    f"{fmt(lbl.x)}\t{fmt(lbl.y)}\t{block}\t{fmt(meshsize)}\t{incircuit}\t"
                    f"{fmt(lbl.magdir)}\t{lbl.group}\t{lbl.turns}\t{ext_default(lbl)}\n"
                )


def add_outer_box(model: FemModel, xmin: float, ymin: float, xmax: float, ymax: float, group: int = 900) -> None:
    pts = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    model.add_polyline(pts, close=True, boundary="A0", group=group)


def add_case1(path: Path, *, variant: str) -> None:
    model = FemModel(comment="Case 1 steel-to-steel return path interruption", depth=20.0)
    model.add_boundary(BoundaryProp("A0"))
    model.add_material(Material("Air", 1.0, 1.0))
    model.add_material(Material("SteelLinear", 4000.0, 4000.0))
    model.add_material(Material("Copper", 1.0, 1.0, sigma=58.0))
    model.add_circuit(Circuit("coil", 1.0))

    add_outer_box(model, -20.0, -20.0, 100.0, 120.0)

    u_core = [(0.0, 0.0), (80.0, 0.0), (80.0, 80.0), (60.0, 80.0),
              (60.0, 20.0), (20.0, 20.0), (20.0, 80.0), (0.0, 80.0)]
    model.add_polyline(u_core, close=True, group=10)
    model.add_segment((20.0, 80.0), (60.0, 80.0), group=11)
    model.add_segment((80.0, 80.0), (80.0, 100.0), group=11)
    model.add_segment((80.0, 100.0), (0.0, 100.0), group=11)
    model.add_segment((0.0, 100.0), (0.0, 80.0), group=11)
    if variant == "truegap":
        gap_box = [(64.0, 49.95), (76.0, 49.95), (76.0, 50.05), (64.0, 50.05)]
        model.add_polyline(gap_box, close=True, group=15)
    else:
        gap_seg = model.add_segment((64.0, 50.0), (76.0, 50.0), group=15)
        if variant == "vgap":
            gap_seg.virtual_gap_enabled = True
            gap_seg.virtual_gap_physical = 0.1
            gap_seg.virtual_gap_mesh = 0.5

    coil = [(24.0, 28.0), (36.0, 28.0), (36.0, 72.0), (24.0, 72.0)]
    model.add_polyline(coil, close=True, group=20)

    model.add_label(BlockLabel(10.0, 10.0, "SteelLinear", meshsize=4.0, group=10))
    model.add_label(BlockLabel(10.0, 90.0, "SteelLinear", meshsize=4.0, group=11))
    model.add_label(BlockLabel(30.0, 50.0, "Copper", meshsize=2.0, circuit="coil", group=20, turns=200))
    model.add_label(BlockLabel(50.0, 50.0, "Air", meshsize=3.0, group=30))
    model.add_label(BlockLabel(-10.0, -10.0, "Air", meshsize=6.0, group=31))
    if variant == "truegap":
        model.add_label(BlockLabel(70.0, 50.0, "Air", meshsize=0.08, group=15))
    model.write(path)


def add_case2(path: Path, *, enable_virtual_gap: bool) -> None:
    model = FemModel(comment="Case 2 pole shoe surface layer / interface perturbation", depth=20.0)
    model.add_boundary(BoundaryProp("A0"))
    model.add_material(Material("Air", 1.0, 1.0))
    model.add_material(Material("SteelLinear", 4000.0, 4000.0))
    model.add_material(Material("Copper", 1.0, 1.0, sigma=58.0))
    model.add_circuit(Circuit("coil", 1.0))

    add_outer_box(model, -20.0, -20.0, 120.0, 110.0)

    base_and_post = [(0.0, 0.0), (100.0, 0.0), (100.0, 20.0), (60.0, 20.0),
                     (60.0, 60.0), (40.0, 60.0), (40.0, 20.0), (0.0, 20.0)]
    model.add_polyline(base_and_post, close=True, group=10)
    model.add_segment((30.0, 60.0), (40.0, 60.0), group=11)
    model.add_segment((60.0, 60.0), (70.0, 60.0), group=11)
    model.add_segment((70.0, 60.0), (70.0, 70.0), group=11)
    model.add_segment((70.0, 70.0), (30.0, 70.0), group=11)
    model.add_segment((30.0, 70.0), (30.0, 60.0), group=11)
    keeper = [(20.0, 70.2), (80.0, 70.2), (80.0, 90.0), (20.0, 90.0)]
    model.add_polyline(keeper, close=True, group=12)
    gap_seg = model.add_segment((32.0, 69.8), (68.0, 69.8), group=16)
    if enable_virtual_gap:
        gap_seg.virtual_gap_enabled = True
        gap_seg.virtual_gap_physical = 0.1
        gap_seg.virtual_gap_mesh = 0.35

    coil_left = [(20.0, 25.0), (34.0, 25.0), (34.0, 55.0), (20.0, 55.0)]
    coil_right = [(66.0, 25.0), (80.0, 25.0), (80.0, 55.0), (66.0, 55.0)]
    model.add_polyline(coil_left, close=True, group=20)
    model.add_polyline(coil_right, close=True, group=21)

    model.add_label(BlockLabel(10.0, 10.0, "SteelLinear", meshsize=4.0, group=10))
    model.add_label(BlockLabel(50.0, 65.0, "SteelLinear", meshsize=1.5, group=11))
    model.add_label(BlockLabel(50.0, 80.0, "SteelLinear", meshsize=3.0, group=12))
    model.add_label(BlockLabel(27.0, 40.0, "Copper", meshsize=2.0, circuit="coil", group=20, turns=100))
    model.add_label(BlockLabel(73.0, 40.0, "Copper", meshsize=2.0, circuit="coil", group=21, turns=-100))
    model.add_label(BlockLabel(-10.0, -10.0, "Air", meshsize=6.0, group=30))
    model.write(path)


def add_case3(path: Path, *, variant: str) -> None:
    model = FemModel(comment="Case 3 permanent magnet to steel interface separation", depth=20.0)
    model.add_boundary(BoundaryProp("A0"))
    model.add_material(Material("Air", 1.0, 1.0))
    model.add_material(Material("SteelLinear", 4000.0, 4000.0))
    model.add_material(Material("PMLinear", 1.05, 1.05, h_c=900000.0))

    add_outer_box(model, -20.0, -20.0, 120.0, 80.0)

    u_frame = [
        (20.0, 0.0),
        (39.5, 0.0),
        (39.5, 38.0),
        (40.5, 38.0),
        (40.5, 0.0),
        (60.0, 0.0),
        (60.0, 50.0),
        (20.0, 50.0),
    ]
    magnet_x0 = 39.5
    magnet = [
        (magnet_x0, 0.0),
        (40.5, 0.0),
        (40.5, 12.0),
        (magnet_x0, 12.0),
    ]

    model.add_polyline(u_frame, close=True, group=10)
    model.add_polyline(magnet, close=True, group=20)
    if variant == "vgap":
        interface_seg = model.add_segment((39.8, 0.5), (39.8, 11.5), group=26)
        interface_seg.virtual_gap_enabled = True
        interface_seg.virtual_gap_physical = 0.1
        interface_seg.virtual_gap_mesh = 0.35

    model.add_label(BlockLabel(26.0, 25.0, "SteelLinear", meshsize=2.5, group=10))
    model.add_label(BlockLabel(40.0, 6.0, "PMLinear", meshsize=0.6, magdir=0.0, group=20))
    model.add_label(BlockLabel(40.0, 22.0, "Air", meshsize=1.0, group=32))
    model.add_label(BlockLabel(-10.0, -10.0, "Air", meshsize=6.0, group=30))
    model.write(path)


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "cases"
    out_dir.mkdir(parents=True, exist_ok=True)
    add_case1(out_dir / "case1_steel_steel_nominal.fem", variant="nominal")
    add_case1(out_dir / "case1_steel_steel_vgap_100um.fem", variant="vgap")
    add_case1(out_dir / "case1_steel_steel_truegap_100um.fem", variant="truegap")
    add_case2(out_dir / "case2_surface_layer_nominal.fem", enable_virtual_gap=False)
    add_case2(out_dir / "case2_surface_layer_vgap_100um.fem", enable_virtual_gap=True)
    add_case3(out_dir / "case3_pm_interface_nominal.fem", variant="nominal")
    add_case3(out_dir / "case3_pm_interface_vgap_100um.fem", variant="vgap")


if __name__ == "__main__":
    main()
