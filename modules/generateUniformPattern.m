function final_result = generateUniformPattern(geometry, target_density, seed, bz_cfg)
% GENERATEUNIFORMPATTERN Applies a stochastic and spatially uniform fibrosis distribution.
%
% Generates a purely random placement of collagen elements. Unlike Perlin-based
% phenotypes, this method does not rely on spatial correlation or biological tensors.
% It guarantees the exact target density based on the Fibrotic Core dimensions.
%
% INPUTS:
%   geometry       - Struct defining valid generation bounds and boundary layers.
%   target_density - Float [0, 1] representing the absolute target collagen fraction.
%   seed           - Integer RNG seed for deterministic generation.
%   bz_cfg         - Struct controlling Gray Zone / Border Zone topological decay rules.
%
% OUTPUTS:
%   final_result   - Struct containing the accumulated outcomes:
%       .presence         - Boolean array of final morphology.
%       .bz_layers        - Topological map inherited from geometry.
%       .fc_density       - Final achieved density within the core.

    rng(seed);

    % Initialize matrices
    target_dims = size(geometry.core_mask);
    presence = false(target_dims);
    
    % Dummy threshold map to maintain compatibility with other modes and visualization
    % 999 represents healthy tissue
    threshold_map = ones(target_dims) * 999; 

    % 1. IDENTIFY TRUE CORE FOR ACCURATE DENSITY CALCULATION
    if isfield(geometry, 'true_core_indices')
        true_core_global_indices = geometry.true_core_indices;
        num_true_core = geometry.num_true_core;
    else
        % Fallback for legacy analytical grids
        true_core_global_indices = geometry.core_indices;
        num_true_core = geometry.num_core_points;
    end

    % 2. Core Generation (Exact density targeting)
    num_target_pixels = round(target_density * num_true_core);
    idx_rand = randperm(num_true_core, num_target_pixels);

    core_pixels = false(num_true_core, 1);
    core_pixels(idx_rand) = true;
    presence(true_core_global_indices) = core_pixels;

    % Viz compatibility: equivalent probability map
    threshold_map(geometry.core_mask) = 1.0 - target_density;

    % 3. Border Zone Generation (Analytical Only)
    if geometry.total_layers > 0
        decay_factor = 3.0; 
        
        for layer = 1:geometry.total_layers
            progression = layer / geometry.total_layers;
            layer_density = target_density * exp(-decay_factor * progression);

            mask_layer = (geometry.bz_map == layer);
            num_layer_points = sum(mask_layer(:));

            if num_layer_points > 0
                num_layer_target = round(layer_density * num_layer_points);
                layer_idx = find(mask_layer);
                idx_rand_layer = randperm(num_layer_points, num_layer_target);

                layer_pixels = false(num_layer_points, 1);
                layer_pixels(idx_rand_layer) = true;

                presence(layer_idx) = layer_pixels;
            end
            threshold_map(mask_layer) = 1.0 - layer_density;
        end
    end

    % 4. Package Result
    final_result.presence = presence;
    final_result.bz_layers = geometry.bz_map;
    final_result.fc_density = sum(presence(true_core_global_indices)) / num_true_core;
end