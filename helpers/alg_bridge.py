"""
alg_bridge.py
Bridge script to handle ALG (CSV-based) I/O between Octave and Python.
"""
import sys
import numpy as np
import scipy.io as sio
import pandas as pd

def write_alg(original_alg, mat_with_presence, out_alg):
    print(f"      |-> Python - Loading framework output data")
    mat_data = sio.loadmat(mat_with_presence)
    presence = mat_data['presence'].flatten()
    
    print(f"      |-> Python - Injecting presence array into ALG structure")
    with open(original_alg, 'r') as fin, open(out_alg, 'w') as fout:
        i = 0
        for line in fin:
            # Skip empty lines to avoid misalignment
            if not line.strip():
                fout.write(line)
                continue
            
            parts = line.split(',')
            if len(parts) >= 7:
                has_collagen = str(int(presence[i]))
                clean_line = line.rstrip('\r\n')
                new_line = f"{clean_line},{has_collagen}\n"
                fout.write(new_line)
                i += 1
            else:
                fout.write(line)
                
    print(f"      |-> [✓] Exported modified mesh to: {out_alg}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(1)
        
    mode = sys.argv[1]
    if mode == 'read':
        print("      |-> [x] Mode not implemented: 'read' is not supported in this script.")
    elif mode == 'write':
        write_alg(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("      |-> [x] Mode not recognized. Use 'read' or 'write'.")