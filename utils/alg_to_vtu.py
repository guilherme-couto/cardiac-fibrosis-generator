"""
Utility script to convert custom ALG meshes into VTU format.
Dynamically handles 2D and 3D column layouts and extracts post-generation fibrosis data.
Supports both Point Cloud (default) and Explicit Voxel (FVM) representations.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import pyvista as pv

def convert_alg_to_vtu(input_csv, output_vtu, is_post_gen, as_cells):
    input_path = Path(input_csv)
    if not input_path.exists():
        print(f"[x] Input file not found: {input_path.resolve()}")
        sys.exit(1)

    if output_vtu is None:
        output_vtu = input_path.with_suffix('.vtu')
    else:
        output_vtu = Path(output_vtu)

    print(f"\n--- Utility: ALG to VTU Converter ---")
    print(f"[i] Reading ALG file: {input_path.name}")
    
    df = pd.read_csv(input_path, sep=',', engine='c', header=None)
    num_cols = df.shape[1]
    
    print(f"[i] Detected {num_cols} columns in the dataset.")
    
    # Extract Spatial Coordinates (Cols 1, 2, 3)
    points = df.iloc[:, 0:3].to_numpy(dtype=np.float32)
    
    # =========================================================================
    # VISUALIZATION MODE SELECTION
    # =========================================================================
    # Option A: Point Data (Default)
    # PROS: Extremely lightweight file size, very fast to load. Perfect for 
    #       rendering biological vectors (using the Glyph filter for arrows).
    # CONS: Renders elements as floating dimensionless points. Cannot visualize 
    #       the solid anatomy as a continuous wall without external filtering.
    #
    # Option B: Cell Data (--as-cells)
    # PROS: Reconstructs the exact Finite Volume Method (FVM) control volumes 
    #       using the dx, dy, dz half-sizes from the ALG file. Renders a true 
    #       continuous solid block anatomy.
    # CONS: File sizes are significantly larger. Generates 8 vertices per element, 
    #       requiring more RAM when opening in ParaView.
    # =========================================================================

    if as_cells:
        print("[i] Mode: CELL DATA. Reconstructing explicit FVM voxels")
        
        # Extract half-dimensions (Cols 4, 5, 6)
        hx = df.iloc[:, 3].to_numpy(dtype=np.float32)
        hy = df.iloc[:, 4].to_numpy(dtype=np.float32)
        hz = df.iloc[:, 5].to_numpy(dtype=np.float32)
        
        x, y, z = points[:, 0], points[:, 1], points[:, 2]

        # Calculate the 8 corners of the hexahedron (voxel) simultaneously
        v0 = np.column_stack((x-hx, y-hy, z-hz))
        v1 = np.column_stack((x+hx, y-hy, z-hz))
        v2 = np.column_stack((x+hx, y+hy, z-hz))
        v3 = np.column_stack((x-hx, y+hy, z-hz))
        v4 = np.column_stack((x-hx, y-hy, z+hz))
        v5 = np.column_stack((x+hx, y-hy, z+hz))
        v6 = np.column_stack((x+hx, y+hy, z+hz))
        v7 = np.column_stack((x-hx, y+hy, z+hz))

        # Vectorized vertex array assembly
        all_verts = np.empty((len(x)*8, 3), dtype=np.float32)
        all_verts[0::8] = v0
        all_verts[1::8] = v1
        all_verts[2::8] = v2
        all_verts[3::8] = v3
        all_verts[4::8] = v4
        all_verts[5::8] = v5
        all_verts[6::8] = v6
        all_verts[7::8] = v7

        # Vectorized cell array assembly for vtkHexahedron
        cells = np.empty((len(x), 9), dtype=int)
        cells[:, 0] = 8 # Number of points per cell
        cells[:, 1] = np.arange(0, len(x)*8, 8)
        for c_idx in range(2, 9):
            cells[:, c_idx] = cells[:, 1] + (c_idx - 1)

        # 12 is the VTK enum for VTK_HEXAHEDRON
        celltypes = np.full(len(x), 12, dtype=np.uint8) 
        
        cloud = pv.UnstructuredGrid(cells.ravel(), celltypes, all_verts)
        target_data = cloud.cell_data
    else:
        print("[i] Mode: POINT DATA. Mapping elements as spatial points")
        cloud = pv.PolyData(points)
        target_data = cloud.point_data

    # 2. Extract Tissue Tag (Col 7)
    if num_cols >= 7:
        tissue_tags = df.iloc[:, 6].to_numpy(dtype=np.int32)
        target_data['Tissue_Tag'] = tissue_tags

    # 3. Extract Longitudinal Fiber (Cols 9, 10, 11)
    if num_cols >= 11:
        fibers = df.iloc[:, 8:11].to_numpy(dtype=np.float32)
        target_data['Fiber_f0'] = fibers
        print("[✓] Extracted longitudinal fiber vectors.")

    # 4. Extract Transverse Sheet and Normal (Cols 12 to 17)
    is_3D_mesh = (num_cols >= 17 and not is_post_gen) or (num_cols >= 18 and is_post_gen)
    if is_3D_mesh:
        sheets = df.iloc[:, 11:14].to_numpy(dtype=np.float32)
        normals = df.iloc[:, 14:17].to_numpy(dtype=np.float32)
        target_data['Sheet_s0'] = sheets
        target_data['Normal_n0'] = normals
        print("[✓] Extracted 3D orthonormal basis (sheet and normal vectors).")
        
    # 5. Extract Laplacian Weights (Second to last column if post-gen, else last column)
    if is_post_gen:
        laplacian_weights = df.iloc[:, -2].to_numpy(dtype=np.float32)
        target_data['Laplacian_Weight'] = laplacian_weights
    else:
        laplacian_weights = df.iloc[:, -1].to_numpy(dtype=np.float32)
        target_data['Laplacian_Weight'] = laplacian_weights
    print("[✓] Extracted Laplacian weights.")

    # 6. Extract Fibrosis Presence (Last Column)
    if is_post_gen:
        collagen = df.iloc[:, -1].to_numpy(dtype=np.int32)
        target_data['Collagen_Presence'] = collagen
        print("[✓] Extracted post-generation collagen presence array.")

    print(f"[i] Saving binary VTU to disk")
    cloud.save(str(output_vtu))
    
    print(f"[✓] Conversion successful! File saved to: {output_vtu.resolve()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts ALG meshes to VTU format.")
    parser.add_argument("input_alg", type=str, help="Path to the input .alg file")
    parser.add_argument("--output", type=str, default=None, help="Optional specific output path for the .vtu file")
    parser.add_argument("--post-gen", action="store_true", help="Indicates the mesh has an appended collagen column.")
    parser.add_argument("--as-cells", action="store_true", help="Reconstructs FVM control volumes (voxels) instead of point clouds.")
    
    args = parser.parse_args()
    
    convert_alg_to_vtu(args.input_alg, args.output, args.post_gen, args.as_cells)