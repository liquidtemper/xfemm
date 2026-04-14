## Virtual Gap Magnetostatics Demos

This folder contains three simple planar magnetostatic benchmark cases for the line-segment Virtual Gap feature.

Cases:

- `case1_steel_steel_*`
  - U-core style electromagnet with a closing steel armature.
  - The virtual gap is placed on a short internal steel-to-steel separation plane in the return path, not on the main pole face.
  - Main metrics:
    - coil flux linkage
    - inductance `L = flux_linkage / current`
    - average normal `B` through the right-leg steel cross section
  - Expected change:
    - the `100 um` gap reduces flux linkage and inductance
    - average `B` in the interrupted return leg decreases
  - Optional true-gap reference:
    - `case1_steel_steel_truegap_100um.fem`
    - implemented as a real `0.1 mm` air slot on the same plane

- `case2_surface_layer_*`
  - Center-post electromagnet with a pole shoe facing a keeper across a `0.5 mm` air gap.
  - The virtual gap is placed at the post-to-pole-shoe interface, acting like a thin nonmagnetic separation layer near the air-facing surface.
  - Main metrics:
    - average normal `B` across the main air-gap center line
    - peak `|B|` in a sampled region around the air gap
    - coil flux linkage / inductance
  - Expected change:
    - the `100 um` layer lowers local permeance
    - average and peak air-gap field decrease

- `case3_pm_interface_*`
  - Permanent magnet closing the bottom opening of a steel U-frame return path.
  - The virtual gap is placed just inside the magnet-side interface region to represent a thin nonmagnetic separation near the magnet-to-steel bond line.
  - Main metrics:
    - average normal `B` through the top steel bridge
    - local `|B|` probes in the steel return path and air cavity
  - Expected change:
    - the interface gap reduces flux driven around the steel return path
    - bridge `B` and local field probes decrease

Dimensions and materials:

- Units: `millimeters`
- Depth: `20 mm`
- Steel: linear `mu_r = 4000`
- Coil current: `1 A`
- Case 1 turns: `200`
- Case 2 turns: `+100 / -100`
- Magnet: linear PM with `mu_r = 1.05`, `Hc = 900 kA/m`
- Virtual-gap physical thickness: `0.1 mm`
- Virtual-gap mesh ribbon thickness:
  - case 1: `0.5 mm`
  - case 2: `0.35 mm`
  - case 3: `0.35 mm`

Automation:

- `generate_virtual_gap_cases.py`
  - regenerates the checked-in `.fem` models
- `run_virtual_gap_benchmarks.py`
  - solves all case variants with `fmesher` + `fsolver`
  - extracts metrics and field samples with `fpproc-vgap`
  - extracts scalar metrics
  - writes a CSV summary
  - samples field magnitude on a grid and writes side-by-side PNG comparisons
  - includes a virtual-gap vs true-gap percent error comparison for case 1

Reference comparison:

- `case1_steel_steel_truegap_100um.fem` is the explicit true-geometry reference model.
- The benchmark summary reports the virtual-gap vs true-gap inductance error for case 1.
