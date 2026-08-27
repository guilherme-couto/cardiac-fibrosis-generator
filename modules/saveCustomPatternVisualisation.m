function saveCustomPatternVisualisation(mesh, output_name, run_cfg)
% SAVECUSTOMPATTERNVISUALISATION Delegates 2D custom mesh rendering to Python.
%
% Triggered when save_figure = true for patient-specific topologies.
% Currently supports 2D .alg and .vtu meshes via scatter plotting.
%
% INPUTS:
%   mesh        - Struct containing framework dimensions (is3D flag).
%   output_name - String containing the full path to the generated output mesh.
%   run_cfg     - Struct providing fibrosis types, target densities, and seeds.

    if mesh.is3D
        fprintf('      |-> [i] Auto-plot is skipped for 3D custom meshes (Use ParaView/alg_to_vtu.py instead).\n');
        return;
    end

    % Extract extension to pass the exact output file
    if strcmpi(mesh.type, 'alg')
        final_mesh_path = [output_name, '.alg'];
    else
        final_mesh_path = [output_name, '.vtu'];
    end

    fibrosis_type = run_cfg.fibrosis_type;
    density = run_cfg.target_density;
    seed = run_cfg.seed;

    title_str = sprintf('Pattern: %s | Density: %.2f | Seed: %d', fibrosis_type, density, seed);
    
    fprintf('      |-> Delegating Custom 2D visual render to Python Matplotlib...\n');
    bridge_script = fullfile(pwd, 'helpers', 'save_custom_fibrosis_image.py');
    
    % Build CLI command
    cmd = sprintf("python3 '%s' '%s' '%s'", bridge_script, final_mesh_path, title_str);
    
    [status, cmd_output] = system(cmd);
    
    if status ~= 0
        fprintf('      |-> [x] Python plotting failed:\n%s\n', cmd_output);
    else
        fprintf('%s', cmd_output); % Prints the success message from Python
    end
end