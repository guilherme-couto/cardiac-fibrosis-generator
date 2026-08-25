function writeALGFromPatternData(pattern_data, filename)
% WRITEALGFROMPATTERNDATA Exports analytical fibrosis patterns to legacy ALG format.
%
% Reconstructs the explicit Cartesian coordinates from the generated analytical 
% spatial arrays and formats them into the finite-volume CSV standard required 
% by MonoAlg3D. Converts internal centimeters back to micrometers for the solver.
%
% INPUTS:
%   pattern_data - Struct encompassing the entire generation result:
%       .mesh      - Base grid information (Nx, Ny, Nz, dx).
%       .presence  - Logical array/matrix of the fibrosis placement.
%       .bz_layers - Topological layer indexing.
%       .params    - Perlin parameters to extract biological fiber angles.
%   filename     - String for the desired output file name (excluding .alg).
%

    Nx = pattern_data.mesh.Nx;
    Ny = pattern_data.mesh.Ny;
    Nz = pattern_data.mesh.Nz;
    is3D = pattern_data.mesh.is3D;
    
    dx_um = pattern_data.mesh.dx * 1e4; % Conversion from cm to micrometers

    fiber_x = cos(pattern_data.params.orientation);
    fiber_y = sin(pattern_data.params.orientation);
    fiber_z = 0;
    
    if is3D
        fiber_x = cos(pattern_data.params.phi) * cos(pattern_data.params.theta);
        fiber_y = -sin(pattern_data.params.phi) * cos(pattern_data.params.theta);
        fiber_z = sin(pattern_data.params.theta);
    end

    num_bz_layers = max(pattern_data.bz_layers(:)) + 1; 

    fprintf('  [i] Exporting analytical array to ALG: %s.alg...\n', filename);
    fid = fopen([filename, '.alg'], 'wt');
    if fid == -1
        error('  [x] Could not create output file. Check directory permissions.'); 
    end

    % Z -> Y -> X ordered iteration
    for k = 1:Nz
        z_center = (k - 0.5) * dx_um;

        for j = 1:Ny
            y_center = (j - 0.5) * dx_um; 

            for i = 1:Nx
                x_center = (i - 0.5) * dx_um;
                
                if is3D
                    fibro_tag = pattern_data.presence(j, i, k); 
                    bz_layer  = pattern_data.bz_layers(j, i, k);
                else
                    fibro_tag = pattern_data.presence(j, i); 
                    bz_layer  = pattern_data.bz_layers(j, i);
                end

                if bz_layer == -1
                    bz_multiplier = 1; 
                else
                    bz_multiplier = bz_layer / num_bz_layers;
                end

                % Export Format: [X, Y, Z, half_dx, half_dy, half_dz, tag, fx, fy, fz, bz_layer, bz_mult]
                fprintf(fid, '%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%d,%.5g,%.5g,%.5g,%d,%.5g\n', ...
                    x_center, y_center, z_center, dx_um*0.5, dx_um*0.5, dx_um*0.5, ...
                    fibro_tag, fiber_x, fiber_y, fiber_z, bz_layer, bz_multiplier);
            end
        end
    end
    fclose(fid);
    fprintf('  [✓] ALG export successful.\n');
end