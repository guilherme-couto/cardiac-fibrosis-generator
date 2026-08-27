"""
Utility script to visualize anatomical tags (Healthy, Core, Gray Zone) 
of a 2D mesh prior to running the fibrosis generator.
"""
import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

def main():
    parser = argparse.ArgumentParser(description="Sanity check for Custom Mesh anatomical tags.")
    parser.add_argument("mesh_file", type=str, help="Path to the .alg or .vtu mesh file")
    args = parser.parse_args()

    filepath = Path(args.mesh_file)
    
    if not filepath.exists():
        print(f"[x] File not found: {filepath.resolve()}")
        sys.exit(1)

    print(f"\n--- Util: Mesh Tag Inspector ---")
    print(f"[i] Reading file: {filepath.name}")

    ext = filepath.suffix.lower()

    if ext == '.vtu':
        grid = pv.read(str(filepath))
        
        if 'material' in grid.cell_data:
            tags = grid.cell_data['material']
        elif 'tecido' in grid.cell_data:
            tags = grid.cell_data['tecido']
        else:
            print("[x] 'material' or 'tecido' array not found in VTU CellData.")
            sys.exit(1)
            
        centers = grid.cell_centers().points
        x, y = centers[:, 0], centers[:, 1]
        print(f"[✓] VTU mesh loaded with {grid.n_cells} elements.")

    elif ext == '.alg':
        df = pd.read_csv(filepath, sep=',', engine='c', header=None)
        x = df.iloc[:, 0].values
        y = df.iloc[:, 1].values
        tags = df.iloc[:, 6].values.astype(int)
        print(f"[✓] ALG mesh loaded with {len(tags)} elements.")
        
    else:
        print(f"[x] Unsupported file extension ({ext}). Use .alg or .vtu.")
        sys.exit(1)

    # --- Plotting Configuration ---
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.spines.top': False,
        'axes.spines.right': False
    })

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    
    # Standard Colors: 0 (Blue), 1 (Orange), 2 (Gray)
    cmap = ListedColormap(['#1A80CC', '#E6801A', '#808080'])
    
    ax.scatter(x, y, c=tags, cmap=cmap, s=5, alpha=0.9, edgecolors='none')

    legend_elements = [
        Patch(facecolor='#1A80CC', label='0: Healthy Tissue'),
        Patch(facecolor='#E6801A', label='1: Fibrotic Core'),
        Patch(facecolor='#808080', label='2: Gray Zone')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, frameon=True)

    ax.set_title(f'Anatomical Tags: {filepath.name}', fontsize=16, pad=15)
    ax.set_xlabel('X Coordinate', fontsize=14)
    ax.set_ylabel('Y Coordinate', fontsize=14)
    ax.axis('equal')

    out_file = f'sanity_check_{filepath.stem}.png'
    plt.savefig(out_file, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"[✓] Visualization saved as: {out_file}\n")

if __name__ == '__main__':
    main()