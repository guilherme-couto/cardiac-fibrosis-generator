"""
vtu_bridge.py
Bridge script to handle VTU I/O via PyVista.
"""
import sys
import pyvista as pv
import scipy.io
import numpy as np

def read_vtu(vtu_file, mat_file):
    print(f"  [i] PyVista - Parsing VTU mesh: {vtu_file}...")
    grid = pv.read(vtu_file)
    
    if 'material' in grid.cell_data:
        tags = grid.cell_data['material']
    elif 'tecido' in grid.cell_data:
        tags = grid.cell_data['tecido']
    else:
        raise ValueError("  [x] Could not find 'material' or 'tecido' arrays in VTU CellData.")
        
    centers = grid.cell_centers().points
    scipy.io.savemat(mat_file, {'points': centers, 'tags': tags})
    print(f"  [✓] Extracted {len(tags)} element centroids.")

def write_vtu(original_vtu, mat_file, output_vtu):
    print(f"  [i] PyVista - Injecting data into VTU geometry...")
    grid = pv.read(original_vtu)
    
    mat_data = scipy.io.loadmat(mat_file)
    new_presence = mat_data['presence'].flatten()
    
    if len(new_presence) != grid.n_cells:
        raise ValueError(f"  [x] Array size mismatch. Presence: {len(new_presence)}, VTU Cells: {grid.n_cells}")
    
    grid.cell_data['material'] = new_presence.astype(np.int32)
    grid.save(output_vtu)
    print(f"  [✓] Exported modified mesh to: {output_vtu}")

if __name__ == "__main__":
    action = sys.argv[1]
    if action == "read":
        read_vtu(sys.argv[2], sys.argv[3])
    elif action == "write":
        write_vtu(sys.argv[2], sys.argv[3], sys.argv[4])