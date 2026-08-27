"""
save_custom_fibrosis_image.py
Matplotlib and PyVista bridge to render 2D Custom Patient Mesh outputs.
Triggered automatically by the Octave framework if save_figure = true.
"""

import sys
from pathlib import Path
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

def plot_custom_mesh(mesh_file, title_str):
    filepath = Path(mesh_file)
    ext = filepath.suffix.lower()

    if ext == '.alg':
        df = pd.read_csv(filepath, sep=',', engine='c', header=None, comment='#')
        x = df.iloc[:, 0].values
        y = df.iloc[:, 1].values
        # The framework appends collagen presence as the very last column
        presence = df.iloc[:, -1].values.astype(int)
    else:
        print(f"      |-> [x] Unsupported extension for auto-plotting: {ext}")
        return

    # --- Plotting Configuration ---
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.spines.top': False,
        'axes.spines.right': False
    })

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    
    color_healthy = '#1A80CC'
    color_fibrosis = '#E6801A'
    cmap = ListedColormap([color_healthy, color_fibrosis])
    
    ax.scatter(x, y, c=presence, cmap=cmap, s=5, alpha=0.9, edgecolors='none')

    legend_elements = [
        Patch(facecolor=color_healthy, label='0: Healthy Tissue'),
        Patch(facecolor=color_fibrosis, label='1: Generated Fibrosis')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, frameon=True)

    ax.set_title(title_str, fontsize=16, pad=15)
    ax.set_xlabel('X Coordinate', fontsize=14)
    ax.set_ylabel('Y Coordinate', fontsize=14)
    ax.axis('equal')

    # Save next to the mesh file
    out_file = filepath.with_name(f"{filepath.stem}.png")
    plt.savefig(out_file, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"      |-> [✓] Custom 2D Image generated: {out_file.name}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("      |-> [x] Usage: python save_custom_fibrosis_image.py <mesh_file> <title>")
        sys.exit(1)
        
    mesh_path = sys.argv[1]
    title = sys.argv[2]
    
    plot_custom_mesh(mesh_path, title)