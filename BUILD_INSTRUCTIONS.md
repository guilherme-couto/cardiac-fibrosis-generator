# WSL Build Instructions and Troubleshooting
**Target System:** Windows Subsystem for Linux (WSL)  
**Core Stack:** Octave 8.4 (Calculation), C++ (MEX Functions), Python (Visualization)  
**Environment Manager:** Conda

## 1. Compilation of C++ Sources
The project contains C++ files (`Octave2D.cpp`, `Octave3D.cpp`) that must be compiled into MEX files callable by Octave. 

**Instructions:**
1. Activate the Conda environment.
2. Navigate to the `helpers/` directory containing the `.cpp` files.
3. Run the following commands:
   ```bash
   mkoctfile --mex Octave2D.cpp
   mkoctfile --mex Octave3D.cpp
   ```

**Expected Output:**
This will generate files with extensions '.mex' (or '.mexa64'). These binaries are automatically detected by Octave as functions.


## 2. KNOWN ISSUE: OCTAVE TERMINAL FREEZE ON WSL
**Description:**
When running Octave via Conda on WSL, the standard 'octave' command may result in a frozen terminal. The prompt accepts no input, and text is not echoed back to the screen.

**Cause:**
This is a conflict between the Conda-provided 'readline' library and the system/WSL interface, often triggered by the GUI (Qt) initialization or terminal capability detection.

**Solution:**
Do not run the standard 'octave' command. Instead, use the command-line interface (CLI) executable with specific flags to disable line editing features.

**Command to run Octave interactively:**
   ```bash
   octave-cli --no-line-editing
   ```

**Command to run scripts directly (Recommended):**
   ```bash
   octave-cli script_name.m
   ```


## 3. VERIFICATION
To verify the build was successful:
   **1. Open the terminal:**
   ```bash
   octave-cli --no-line-editing
   ```
   **2. Check for MEX files:**
   ```bash
   ls *.mex
   ```
   **3. Attempt to call the function helper (even if undocumented):**
   ```bash
   help Octave2D
   ```

If the error returns "'Octave2D' is not documented" rather than "undefined," the binary has been successfully loaded.


## 4. THE PYTHON FOR IMAGES
**Problem:**
The Octave binary provided by Conda often has Image I/O disabled at compile time, causing "support unavailable" errors.

**Solution:**
Visualization is handled by an external Python script (`helpers/save_fibrosis_image.py`).
   1. Octave calculates the fibrosis pattern.
   2. Octave saves the raw data matrix to a temporary .mat file.
   3. Octave triggers the Python script via system command.
   4. Python (Matplotlib) reads the .mat file and saves the final PNG.

**Requirements:**
   1. Ensure `save_fibrosis_image.py` is in the `helpers/` directory.
   2. Ensure `matplotlib`, `scipy`, `numpy` and `pyvista` are installed in the Conda environment.


## 5. INSTALLING OCTAVE PACKAGES
If the project requires standard Octave toolboxes (e.g., `image` package) for calculation logic:

**Prerequisite:** Ensure `graphicsmagick` is installed in Conda.

**Installation:**
   Open Octave and run:
   ```bash
   pkg install -forge image
   pkg load image
   ```

**Note:** This requires Octave 8.4.0. Newer versions (10.x) may fail to compile due to GCC incompatibilities


## 6. TROUBLESHOOTING
**Error:** 'timespec_get' has not been declared during package installation.
   **Fix:** You are likely using Octave 10 with GCC 15. Downgrade to Octave 8.4.0 using Conda.

**Error:** imwrite: support for Image IO was unavailable.
   **Fix:** Do not use imwrite in Octave code. Use the Python bridge method described in Section 4.