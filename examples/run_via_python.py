# EXAMPLE: PYTHON ORCHESTRATOR
# Demonstrates how to call the Cardiac Fibrosis Framework from Python.
# This script is the direct Python equivalent of the Octave example scripts.

import subprocess
from pathlib import Path

print("Cardiac Fibrosis Framework Example: Python Orchestration\n")

# Helper function to safely format Python variables into Octave syntax
def fmt(arg):
    if isinstance(arg, str): return f"'{arg}'"
    if isinstance(arg, list): return "[" + ", ".join(map(str, arg)) + "]"
    if arg is None or (isinstance(arg, list) and len(arg) == 0): return "[]"
    if isinstance(arg, bool): return "true" if arg else "false"
    return str(arg)

# Configure paths (Relative to the framework root)
filepath_2D = '../path/to/your/custom_mesh.alg'  # Replace with your actual .alg file path

# --- SCENARIO 1: Analytical Grid ---
scn1 = {
    'id': 5,
    'type': 'interstitial',
    'dim_mode': '2D',
    'shape': 'ellipse',
    'angle': 60,
    'density': 0.10,
    'domain': [0.01, 2.0, 4.0],  # [dx, Lx, Ly]
    'core': [1.5, 3.0],          # [width, height]
    'desc': '2D Elliptical Domain - Interstitial'
}

# --- SCENARIO 2: Custom Patient Mesh ---
scn2 = {
    'id': 6,
    'type': 'diffuse',
    'dim_mode': 'CUSTOM',
    'shape': [],
    'angle': [],
    'density': 0.40,
    'domain': filepath_2D,       # String path
    'core': [],
    'desc': 'Patient (2D Slice) - Diffuse'
}

test_scenarios = [scn1, scn2]
seed = 2026
save_mesh = True
save_figure = True

# Navigate to the root directory of the framework (one level up from /examples)
root_dir = Path(__file__).resolve().parent.parent

# --- EXECUTION LOOP ---
for i, s in enumerate(test_scenarios, 1):
    
    fname = f"id{s['id']}_{s['dim_mode']}_{s['type']}_{seed}"
                
    print(f">>> Running Scenario {i}: {s['desc']}")
    
    # Define output path
    output_path = Path('examples/outputs') / fname
    
    # Build the exact Octave function call
    octave_cmd = (
        f"run_fibrosis_generator({fmt(s['type'])}, {fmt(s['density'])}, {fmt(seed)}, "
        f"{fmt(s['angle'])}, {fmt(s['dim_mode'])}, {fmt(s['domain'])}, {fmt(s['shape'])}, "
        f"{fmt(s['core'])}, {fmt(str(output_path))}, {fmt(save_mesh)}, {fmt(save_figure)});"
    )
    
    # Mount the subprocess command
    cmd = ["octave-cli", "-W", "-q", "--no-gui", "--eval", octave_cmd]
    
    try:
        # Execute Octave inside the framework's root directory
        subprocess.run(cmd, check=True, cwd=str(root_dir))
    except subprocess.CalledProcessError as e:
        print("[x] Error executing Octave via Python.\n")