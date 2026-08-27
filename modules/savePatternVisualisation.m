function savePatternVisualisation(output_data, final_result, run_cfg)
% SAVEPATTERNVISUALISATION Delegates rapid graphical rendering to Python.
%
% Triggers an external headless Python script (save_fibrosis_image.py) to parse 
% the structural matrix and generate high-quality standardized PNG visualizations 
% of the generated analytical geometries, bypassing Octave's native plotting limits.
%
% INPUTS:
%   output_data  - Struct containing the overarching parameters and 3D flag.
%   final_result - Struct containing the targeted .presence array matrix.
%   run_cfg      - Struct providing output filenames, fibrosis types, and densities.

    presence = final_result.presence;
    
    output_name = run_cfg.output_name;
    fibrosis_type = run_cfg.fibrosis_type;
    density = run_cfg.target_density;
    seed = run_cfg.seed;

    temp_mat_file = [output_name, '_temp.mat'];
    save("-v7", temp_mat_file, "presence");
    
    if isfield(output_data, 'is3D') && output_data.is3D
        dim_tag = '3D';
    else
        dim_tag = '2D';
    end
    title_str = sprintf('%s (%s) | Density: %.2f | Seed: %d', fibrosis_type, dim_tag, density, seed);
    
    fprintf('      |-> Delegating %s visual render to Python Matplotlib\n', dim_tag);
    bridge_script = fullfile(pwd, 'helpers', 'save_fibrosis_image.py');
    cmd = sprintf("python3 '%s' '%s' '%s' '%s'", bridge_script, temp_mat_file, output_name, title_str);
    
    [status, cmd_output] = system(cmd);
    
    if status ~= 0
        fprintf('      |-> [x] Python plotting failed:\n%s\n', cmd_output);
    else
        fprintf('      |-> [✓] Image generated.\n');
    end
end