"""
Pre-processing utility to calculate spatial Gray Zone weights for 2D meshes.
Solves a steady-state Laplacian PDE to map the continuous transition 
between the Fibrotic Core (0.0) and Healthy Tissue (1.0).
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def compute_grayzone_laplacian(input_csv, output_csv, img_output):
    input_path = Path(input_csv)
    if not input_path.exists():
        print(f"[x] Input file not found: {input_path.resolve()}")
        sys.exit(1)

    print(f"\n--- Utility: 2D Laplacian Weight Solver ---")
    print(f"[i] Reading mesh: {input_path.name}")
    
    df = pd.read_csv(input_path, sep=',', header=None)
    
    points = df.iloc[:, 0:3].values
    tags = df.iloc[:, 6].values
    
    idx_health = np.where(tags == 0)[0]
    idx_core = np.where(tags == 1)[0]
    idx_gz = np.where(tags == 2)[0]
    
    if len(idx_gz) == 0:
        print("[!] No Gray Zone (Tag 2) found in the mesh.")
        print(f"[i] Copying original file directly to {output_csv}")
        df.to_csv(output_csv, sep=',', header=False, index=False)
        return

    print("[i] Building structured grid indices (von Neumann neighborhood)")
    X = np.round(points[:, 0], 6)
    Y = np.round(points[:, 1], 6)
    
    unique_x = np.unique(X)
    unique_y = np.unique(Y)
    
    x_to_i = {val: i for i, val in enumerate(unique_x)}
    y_to_j = {val: j for j, val in enumerate(unique_y)}
    
    coord_to_idx = {}
    for idx in range(len(points)):
        i = x_to_i[X[idx]]
        j = y_to_j[Y[idx]]
        coord_to_idx[(i, j)] = idx
        
    num_points = len(points)
    A = lil_matrix((num_points, num_points))
    b = np.zeros(num_points)
    
    print("[i] Assembling sparse Laplacian matrix")
    for idx in range(num_points):
        if tags[idx] == 0:   
            # Healthy Tissue -> Dirichlet boundary: W = 1.0
            A[idx, idx] = 1.0
            b[idx] = 1.0
        elif tags[idx] == 1: 
            # Fibrotic Core -> Dirichlet boundary: W = 0.0
            A[idx, idx] = 1.0
            b[idx] = 0.0
        else:                
            # Gray Zone -> Laplacian PDE
            i = x_to_i[X[idx]]
            j = y_to_j[Y[idx]]
            
            neighbors = []
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                ni, nj = i + di, j + dj
                if (ni, nj) in coord_to_idx:
                    neighbors.append(coord_to_idx[(ni, nj)])
                    
            A[idx, idx] = len(neighbors)
            for n_idx in neighbors:
                A[idx, n_idx] = -1.0
                
    print("[i] Solving PDE via scipy.sparse")
    weights = spsolve(A.tocsr(), b)
    
    # Clip to avoid floating point overshoot
    weights = np.clip(weights, 0.0, 1.0)
    
    print(f"[i] Appending weights to ALG and saving")
    df[df.shape[1]] = weights 
    df.to_csv(output_csv, sep=',', header=False, index=False)
    
    print("[i] Generating validation plot")
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.spines.top': False,
        'axes.spines.right': False
    })

    plt.figure(figsize=(10, 8), dpi=300)
    sc = plt.scatter(points[:, 0], points[:, 1], c=weights, cmap='coolwarm', s=10, marker='s')
    cbar = plt.colorbar(sc)
    cbar.set_label('Laplacian Weight (0 = Core, 1 = Healthy)', fontsize=12)
    
    plt.title(f'Grayzone PDE Weights: {input_path.name}', fontsize=16, pad=15)
    plt.xlabel('X Coordinate', fontsize=14)
    plt.ylabel('Y Coordinate', fontsize=14)
    plt.axis('equal')
    
    plt.savefig(img_output, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"[✓] Pre-processing complete!")
    print(f"[✓] Output Mesh: {output_csv}")
    print(f"[✓] Validation Image: {img_output}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes Gray Zone Laplacian weights for 2D meshes.")
    parser.add_argument("input_alg", type=str, help="Path to the input .alg file")
    parser.add_argument("output_alg", type=str, help="Path to save the output .alg file")
    parser.add_argument("--img", type=str, default="laplacian_weights.png", help="Path to save the validation image (default: laplacian_weights.png)")
    
    args = parser.parse_args()
    
    compute_grayzone_laplacian(args.input_alg, args.output_alg, args.img)