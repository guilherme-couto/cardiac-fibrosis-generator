function final_result = runCompositionLoop(mesh, geometry, noise_params, target_density, bz_cfg)
% RUNCOMPOSITIONLOOP Iterative algorithm for sparse procedural accumulation.
%
% Progressively adds fractal noise layers, evaluating spatial occupancy strictly 
% inside the anatomical True Core. The boolean union of these fields forms 
% highly disconnected, patchy morphologies until the target density is reached.
%
% INPUTS:
%   mesh           - Struct containing spatial coordinates and biological tensors.
%   geometry       - Struct defining valid generation bounds and boundary layers.
%   noise_params   - Struct with Perlin variables (seed, fibreness, patchiness).
%   target_density - Float [0, 1] representing the absolute target collagen fraction.
%   bz_cfg         - Struct controlling Gray Zone topological decay rules.
%
% OUTPUTS:
%   final_result   - Struct containing the accumulated outcomes:
%       .presence         - Boolean array of final morphology.
%       .bz_layers        - Topological map inherited from geometry.
%       .fc_density       - Final achieved density within the core.

    target = target_density;
    current_density = 0;
    
    if isfield(mesh, 'tags')
        accumulated_presence = false(mesh.num_points, 1);
    else
        if mesh.is3D
            accumulated_presence = false(mesh.Ny, mesh.Nx, mesh.Nz);
        else
            accumulated_presence = false(mesh.Ny, mesh.Nx);
        end
    end

    step_result.presence = false(size(accumulated_presence));
    step_result.threshold_map = ones(size(geometry.core_mask)) * 999;

    min_step = 1.0 / geometry.num_true_core; 
    tolerance = max(0.005, 1.5 * min_step);
    
    req_density = target;
    iter = 0;
    max_tries = 10;
    tries = 0;
        
    while abs(current_density - target) > tolerance
        iter = iter + 1;
        
        seed_iter = noise_params.seed + iter;
        noise_field = computeNoiseField(mesh, noise_params, seed_iter);
        
        if req_density >= 0.1 - tolerance
             eff_thresh_dens = 0.1; 
        else
             eff_thresh_dens = req_density;
        end
        
        step_result = applyFibrosisThreshold(noise_field, geometry, eff_thresh_dens, bz_cfg);
        temp_presence = accumulated_presence | step_result.presence;
        temp_density = sum(temp_presence(geometry.true_core_indices)) / geometry.num_true_core;
        
        tries = tries + 1;
        
        if temp_density < target + tolerance
            accumulated_presence = temp_presence;
            current_density = temp_density;
            tries = 0;
            
            req_density = abs(current_density - target);
        end
        
        if tries > max_tries
            fprintf('      [i] Composition loop settled at stable state. Density: %.2f%% (Target: %.2f%%)\n', current_density * 100, target * 100);
            break; 
        end
        
        if iter > 200 
            fprintf('      [i] Maximum composition steps (200) reached. Density: %.2f%% (Target: %.2f%%)\n', current_density * 100, target * 100);
            break; 
        end
    end
    
    if abs(current_density - target) <= tolerance
        fprintf('      [i] Composition convergence reached smoothly. Density: %.2f%%\n', current_density * 100);
    end
    
    final_result.presence = accumulated_presence;
    final_result.bz_layers = geometry.bz_map;
    final_result.fc_density = current_density;
end