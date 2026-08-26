# Cardiac Fibrosis Generator

A computational framework for generating synthetic histology-inspired fibrosis patterns (e.g., Compact, Diffuse, Interstitial, Patchy, Uniform) in cardiac models. 

This tool supports both **analytical ideal grids** and **patient-specific geometries**.

---

## Prerequisites & Installation

The framework relies on a hybrid Octave/Python architecture for high-performance generation and visualization.

* **Calculation Engine:** Octave 8.4+ (requires the `image` package: `pkg install -forge image`).
* **Visuals & Large I/O:** Python 3.x (`pyvista`, `pandas`, `scipy`, `matplotlib`).
* **C++ Compiler:** Required to build the core noise functions into MEX binaries.

> **Important:** Before running any scripts, you must compile the C++ MEX binaries. Please read [BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md) for the exact compilation steps and WSL troubleshooting.

---

## Quick Start

The framework is operated via the universal entry point: `run_fibrosis_generator.m`.
To test the installation, run the provided examples from the `examples/` directory:

```bash
cd examples

% Example 1: Generate analytical grids (2D/3D)
octave-cli run_ideal_meshes.m

% Example 2: Patterns into custom patient meshes
octave-cli run_custom_meshes.m

% Example 3: Run the generator via Python (2D/3D/custom)
python3 run_via_python.py
```

---

## Basic Functioning & Density Rule

The algorithm places collagen elements based on a generated continuous noise field. It automatically maps the available space (Core vs. Healthy) and thresholds the noise to achieve the requested volumetric density.

### Crucial Rule: Density Targeting
**The target density requested by the user applies only to the Fibrotic Core.**
Collagen elements that fall into the anatomical Gray Zone (in custom meshes) or Border Zone layers (in analytical meshes) are generated as a decaying byproduct and **do not** count towards the user-requested target density. 

*Example: Requesting a density of `0.4` means exactly 40% of the Fibrotic Core will become collagen, regardless of how much collagen spills into the surrounding border zones.*

---

## Custom Mesh Requirements (.ALG & .VTU)

When using the `CUSTOM` dimension mode, your patient-specific meshes must follow some formatting rules so the framework can correctly extract spatial coordinates, biological tensors, and anatomy.

**Tissue Tagging Convention (Both Formats):**
* `0` = Healthy Tissue
* `1` = Fibrotic Core (Dense target)
* `2` = Gray Zone / Border Zone (Decay target)

### ALG Format (CSV-based Finite Volume)
The framework utilizes a parser that targets specific columns (1-indexed) within the `.alg` file. Unused columns are safely ignored during generation but perfectly preserved in the final output.

* **Columns 1, 2, 3:** X, Y, Z spatial coordinates of centroids.
* **Column 7:** Tissue Tag (0, 1, or 2).
* **Columns 9, 10, 11:** Longitudinal Fiber vector (fx, fy, fz).
* **Columns 12, 13, 14 (3D Only):** Transverse Sheet vector (sx, sy, sz).
* **Columns 15, 16, 17 (3D Only):** Sheet Normal vector (nx, ny, nz).
* **Column 12 (for 2D) or 18 (for 3D):** Pre-calculated Gray Zone Laplacian weights [0.0 to 1.0].

### VTU Format
For `.vtu` grids, spatial centroids are extracted automatically via PyVista. However, the file **must** contain a `CellData` array named exactly `material` or `tecido` holding the integer tissue tags (0, 1, or 2).

---

## Key Customizations

You can tweak the topological behavior of the generated meshes directly inside the module files:

### 1. Gray Zone / Border Zone Decay
Located in `modules/applyFibrosisThreshold.m`.
You can alter how quickly the collagen presence decays outside the true core by changing the `gz_decay_mode`:
* **`linear`**: Standard linear interpolation. Causes a rapid drop in density towards the healthy tissue.
* **`power`**: Uses a fractional power equation ($T^{1-W^k}$). This shields the Gray Zone, preserving a higher density of collagen further away from the core before decaying. You can adjust the intensity via the `power_k` variable.

---

## Related Papers

> Couto, G. M., Campos, J. O., & dos Santos, R. W. (2026). A Framework to Assess the Influence of Different Cardiac Fibrosis Phenotypes and Border Zone on Arrhythmogenesis. In *Simpósio Brasileiro de Computação Aplicada à Saúde (SBCAS)* (pp. 1463-1468). SBC.  
> DOI: [10.5753/sbcas.2026.21656](https://doi.org/10.5753/sbcas.2026.21656)

> Lawson, B. A., Drovandi, C., Burrage, P., Bueno-Orovio, A., Dos Santos, R. W., Rodriguez, B., Mengersen, K., & Burrage, K. (2024). Perlin noise generation of physiologically realistic cardiac fibrosis. *Medical Image Analysis*, 98, 103240.  
> DOI: [10.1016/j.media.2024.103240](https://doi.org/10.1016/j.media.2024.103240)