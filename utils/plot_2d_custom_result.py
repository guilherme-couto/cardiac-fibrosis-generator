"""
Utility script to generate a 2D plot visualization 
of the generated fibrosis pattern on Custom Meshes (.alg).
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

def main():
    parser = argparse.ArgumentParser(description="Plots the generated fibrosis pattern on a 2D Custom ALG mesh.")
    parser.add_argument("mesh_file", type=str, help="Path to the post-generation .alg file")
    parser.add_argument("--dpi", type=int, default=300, help="Output image DPI (default: 300)")
    args = parser.parse_args()

    filepath = Path(args.mesh_file)
    
    if not filepath.exists():
        print(f"[x] File not found: {filepath.resolve()}")
        sys.exit(1)

    print(f"\n--- Utility: 2D Custom Mesh Result Plotter ---")
    print(f"[i] Reading file: {filepath.name}")

    if filepath.suffix.lower() != '.alg':
        print(f"[x] This script currently supports post-generation .alg files only.")
        sys.exit(1)

    # 1. Load the data using Pandas C-engine for speed
    df = pd.read_csv(filepath, sep=',', engine='c', header=None, comment='#')
    
    # 2. Extract Coordinates
    x = df.iloc[:, 0].values
    y = df.iloc[:, 1].values
    
    # 3. Extract Fibrosis Presence (The framework appends it as the very last column)
    presence = df.iloc[:, -1].values.astype(int)
    
    print(f"[✓] Extracted {len(presence)} elements.")
    print(f"     |-> Healthy (0): {sum(presence == 0)}")
    print(f"     |-> Fibrosis (1): {sum(presence == 1)}")

    # --- Plotting Configuration ---
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.spines.top': False,
        'axes.spines.right': False
    })

    fig, ax = plt.subplots(figsize=(8, 8), dpi=args.dpi)
    
    # Standard Fibrosis Colors
    color_healthy = '#1A80CC'
    color_fibrosis = '#E6801A'
    cmap = ListedColormap([color_healthy, color_fibrosis])
    
    # Scatter plot for unstructured points
    ax.scatter(x, y, c=presence, cmap=cmap, s=5, alpha=0.9, edgecolors='none')

    legend_elements = [
        Patch(facecolor=color_healthy, label='0: Healthy Tissue'),
        Patch(facecolor=color_fibrosis, label='1: Fibrosis')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, frameon=True)

    # Use the filename without extension as the title
    ax.set_title(f'Procedural Pattern: {filepath.stem}', fontsize=16, pad=15)
    ax.set_xlabel('X Coordinate', fontsize=14)
    ax.set_ylabel('Y Coordinate', fontsize=14)
    ax.axis('equal')

    out_file = f'render_{filepath.stem}.png'
    plt.savefig(out_file, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"[✓] Visualization saved as: {out_file}\n")

if __name__ == '__main__':
    main()