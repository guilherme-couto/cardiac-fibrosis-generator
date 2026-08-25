function result = applyFibrosisThreshold(noise_field, geometry, target_density, bz_cfg)
% APPLYFIBROSISTHRESHOLD Maps a continuous noise field to a boolean presence array.
%
% Evaluates the generated procedural noise over the defined geometry, using a 
% binary search to find the exact threshold that achieves the requested target 
% density within the True Core. Also computes the Gray Zone (Border Zone) decay.
%
% INPUTS:
%   noise_field    - [N x 1] Float array containing normalized Perlin noise [0, 1].
%   geometry       - Struct defining the core mask, true core indices, and border zones.
%   target_density - Float [0, 1] specifying the desired collagen volume fraction.
%   bz_cfg         - Struct containing Border Zone configurations (metric, layers).
%
% OUTPUTS:
%   result         - Struct containing:
%       .presence      - Boolean array/matrix of finalized fibrosis placement.
%       .bz_layers     - Array/matrix tracking border zone layer depth.
%       .fc_density    - Actual achieved density purely inside the true core.
%       .threshold_map - Spatial map of the explicit threshold values applied.

    noise_field = noise_field(:);
    max_iters = 40;
    tol = 1e-3;

    % 1. Identify True Core
    if isfield(geometry, 'weights') && ~isempty(geometry.weights)
        true_core_local_mask = (geometry.tags == 1);
        num_true_core = sum(true_core_local_mask);
        
        % Fallback if the region is entirely composed of Gray Zone
        if num_true_core == 0
            true_core_local_mask = (geometry.tags == 2);
        end
        true_core_global_indices = geometry.core_indices(true_core_local_mask);
        core_values = noise_field(true_core_global_indices);
    else
        core_values = noise_field(geometry.core_indices);
    end

    % 2. Binary Search for Base Threshold
    [core_threshold, actual_core_density] = findThresholdForDensity(core_values, target_density, tol, max_iters);

    % 3. Build Threshold Map with Gray Zone Decay
    threshold_map = ones(size(noise_field)) * 999; 
    
    if isfield(geometry, 'weights') && ~isempty(geometry.weights)
        
        % --- GRAY ZONE DECAY CONFIGURATION ---
        % 'linear' : Classic linear interpolation (fast drop in density)
        % 'power'  : Power function T^(1-W^k) (preserves GZ density longer)
        gz_decay_mode = 'power'; 
        power_k = 2.0; % Adjust this to control the penalty intensity (k > 1 protects GZ)
        
        W = geometry.weights;
        
        if strcmpi(gz_decay_mode, 'power')
            % Mathematical Shielding: 
            % Base of a fractional power cannot be negative (avoids complex numbers)
            % and cannot be exactly zero (0^0 is NaN).
            T_safe = max(1e-6, min(0.999, core_threshold));
            
            % If threshold is astronomically high (0% density) or low (100% density),
            % we bypass the equation to maintain the absolute block/pass.
            if core_threshold > 1.0
                local_thresholds = ones(size(W)) * 999;
            elseif core_threshold < 0.0
                local_thresholds = ones(size(W)) * -999;
            else
                % The custom power equation: T_layer = T_core^(1 - W^k)
                local_thresholds = T_safe .^ (1.0 - (W .^ power_k));
            end
        else
            % Standard Linear Equation: T_local = T_core + W * (1 - T_core)
            local_thresholds = core_threshold + W .* (1.0 - core_threshold);
        end
        
        threshold_map(geometry.core_indices) = local_thresholds;
    else
        threshold_map(geometry.core_mask) = core_threshold;
    end
    
    % Maintains support for classical analytical Border Zone layers
    if isfield(geometry, 'total_layers') && geometry.total_layers > 0
        decay_factor = 3.0; 
        L = geometry.total_layers;
        for layer = 1:L
            mask_layer = (geometry.bz_map == layer);
            layer_values = noise_field(mask_layer);
            if isempty(layer_values), continue; end
            
            target_layer_density = target_density * exp(-decay_factor * (layer / L));
            layer_thresh = findThresholdForDensity(layer_values, target_layer_density, tol, max_iters);
            threshold_map(mask_layer) = layer_thresh;
        end
    end
    
    % 4. Final Application
    presence = (noise_field >= threshold_map);

    target_dims = size(geometry.core_mask);
    result.presence = reshape(presence, target_dims);
    
    if isfield(geometry, 'bz_map')
        result.bz_layers = geometry.bz_map;
    else
        result.bz_layers = [];
    end
    
    result.fc_density = actual_core_density;
    result.threshold_map = reshape(threshold_map, target_dims);
end

% Helper Function
function [thresh, actual_dens] = findThresholdForDensity(values, target_dens, tol, max_iters)
    if target_dens <= 0
        thresh = max(values) + 0.1; 
        actual_dens = 0;
        return;
    end
    if target_dens >= 1
        thresh = min(values) - 0.1; 
        actual_dens = 1;
        return;
    end

    low = min(values);
    high = max(values);
    
    for i = 1:max_iters
        thresh = (low + high) / 2;
        actual_dens = sum(values >= thresh) / length(values);
        
        if abs(actual_dens - target_dens) < tol
            break;
        end
        if actual_dens < target_dens
            high = thresh; 
        else
            low = thresh;  
        end
    end
end