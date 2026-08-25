function writeCustomMesh(mesh, presence, output_filename)
% WRITECUSTOMMESH Injects procedural arrays into patient-specific mesh files.
%
% Bridges Octave and Python to securely rewrite massive unstructured datasets. 
% Saves the generated boolean presence array as a temporary binary buffer and 
% calls an external Python bridge to append the data as a new 
% column to the original .alg or .vtu file, preserving the anatomical header.
%
% INPUTS:
%   mesh            - Struct containing .type ('alg'/'vtu') and .original_file path.
%   presence        - [N x 1] Boolean or float array representing collagen.
%   output_filename - String specifying the desired output name (without extension).

    fprintf('  [i] Exporting procedural pattern to patient mesh...\n');

    if strcmpi(mesh.type, 'alg')
        
        out_file = [output_filename, '.alg'];
        temp_mat = 'temp_alg_out.mat';
        
        save("-v7", temp_mat, "presence");
        
        bridge_script = fullfile(pwd, 'helpers', 'alg_bridge.py');
        cmd = sprintf("python3 '%s' write '%s' '%s' '%s'", bridge_script, mesh.original_file, temp_mat, out_file);
                      
        [status, cmdout] = system(cmd);
        
        if status ~= 0
            error('  [x] Python ALG bridge failed:\n%s', cmdout);
        end

        delete(temp_mat);
        
    elseif strcmpi(mesh.type, 'vtu')
        out_file = [output_filename, '.vtu'];
        temp_mat = 'temp_vtu_out.mat';
        
        save('-v7', temp_mat, 'presence');
        
        bridge_script = fullfile(pwd, 'helpers', 'vtu_bridge.py');
        cmd = sprintf("python3 '%s' write '%s' '%s' '%s'", bridge_script, mesh.original_file, temp_mat, out_file);
        
        [status, cmd_out] = system(cmd);
        if status ~= 0
            error('  [x] Python VTU bridge failed:\n%s', cmd_out);
        end
        
        delete(temp_mat);
    end
end