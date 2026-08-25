function [combined_noise, noise_comps] = computeNoiseField(mesh, params, seed)
% COMPUTENOISEFIELD Generates a normalized procedural noise field.
%
% Blends multiple structural noise components (base noise, fibrous clefts, 
% and density variations) advected by the local biological orthonormal basis 
% to produce physiologically realistic textures.
%
% INPUTS:
%   mesh   - Struct containing spatial coordinates (.mm_points) and tensors.
%   params - Struct with procedural variables (fibreness, feature_size, etc.).
%   seed   - Integer RNG seed for deterministic generation.
%
% OUTPUTS:
%   combined_noise - [N x 1] Array of the finalized continuous noise [0, 1].
%   noise_comps    - Struct containing the isolated components (O_b, F, O_d).

    % Check Dimensionality
    is3D = mesh.is3D;

    % === COMMON PARAMETER EXTRACTION ===
    fibreness    = params.fibreness;
    patchiness   = params.patchiness;
    feature_size = params.feature_size;
    roughness    = params.roughness;
    patch_size   = params.patch_size;

    % Generate Tables
    if nargin < 3, seed = params.seed; end
    [perm_table, offsets1, offsets2] = generateTables(seed, is3D);

    % Precompute perm tables for structured noise
    perm_table2 = zeros(size(perm_table), 'int32');
    perm_table3 = zeros(size(perm_table), 'int32');
    perm_table4 = zeros(size(perm_table), 'int32');
    
    for k = 1:size(perm_table,1)
        perm_table2(k,:) = perm_table(k, perm_table(k,:)+1);
        perm_table3(k,:) = perm_table(k, perm_table2(k,:)+1);
        if is3D
            perm_table4(k,:) = perm_table(k, perm_table3(k,:)+1);
        end
    end

    % =====================================================================
    % BRANCH: 3D GENERATION (FLOW NOISE UPGRADE + COSINE VEINS)
    % =====================================================================
    if is3D
        % Extract 3D Specific Params
        fibre_sep_y = params.fibre_sep_y;
        fibre_sep_z = params.fibre_sep_z;
        align_y     = params.alignment_y;
        align_z     = params.alignment_z;
        phi         = params.phi;
        theta       = params.theta;

        % Handle Non-Fibrous Case
        if isnan(fibreness)
            fibreness = 0; fibre_sep_y = 1; fibre_sep_z = 1;
        end

        points = mesh.mm_points; % Nx3 Matrix
        P_base = points' / feature_size; % Transpose to 3xN format for MEX engine
        
        if isfield(mesh, 'fibers')
            % Extract and normalize the primary fiber vector
            f0 = mesh.fibers;
            f0 = f0 ./ max(sqrt(sum(f0.^2, 2)), 1e-8);
            
            % Ensure the 3D mesh possesses a complete orthonormal basis (Sheet & Normal)
            if isfield(mesh, 'sheet') && isfield(mesh, 'normal')
                s0 = mesh.sheet;
                n0 = mesh.normal;
                s0 = s0 ./ max(sqrt(sum(s0.^2, 2)), 1e-8);
                n0 = n0 ./ max(sqrt(sum(n0.^2, 2)), 1e-8);
            else
                % Fallback synthesis if the 3D mesh only provides the longitudinal fiber
                tmp = [ones(size(f0,1),1), zeros(size(f0,1),2)];
                s0 = cross(f0, tmp, 2);
                invalid = sum(s0.^2, 2) < 1e-5;
                s0(invalid, :) = cross(f0(invalid,:), [zeros(sum(invalid),1), ones(sum(invalid),1), zeros(sum(invalid),1)], 2);
                s0 = s0 ./ max(sqrt(sum(s0.^2, 2)), 1e-8);
                n0 = cross(f0, s0, 2);
            end
            
            % Apply global rotation offsets (phi/theta)
            R_offset = [ cos(phi) * cos(theta), sin(phi), -cos(phi) * sin(theta);
                        -cos(theta) * sin(phi), cos(phi),  sin(phi) * sin(theta);
                         sin(theta)           , 0       ,  cos(theta) ];
            
            f0_rot = (R_offset * f0')';
            s0_rot = (R_offset * s0')';
            n0_rot = (R_offset * n0')';
            
            % STEP 1: MAIN NOISE FIELD (O_b) - Advected natively by C++
            O_b = Octave3D(P_base, 4, roughness, perm_table, offsets1, f0_rot', s0_rot', n0_rot', align_y, align_z);

            % STEP 2: FIBRE-SELECTING FIELD (F) - The interstitial cleft generator
            centroid = mean(points, 1);
            v = points - centroid;
            y_prime = sum(v .* s0_rot, 2)'; % Transversal axis 1 (1xN)
            z_prime = sum(v .* n0_rot, 2)'; % Transversal axis 2 (1xN)
            
            n_fibres_similarity = 4; wiggle_feature_length = 4; phasefield_strength = 5;
            P_veins = points' / wiggle_feature_length;
            
            phasefield1 = Octave3D(P_veins, 4, 0.5, perm_table2, offsets1, f0_rot', s0_rot', n0_rot', 1.0, 1.0);
            phasefield2 = Octave3D(P_veins, 4, 0.5, perm_table3, offsets2, f0_rot', s0_rot', n0_rot', 1.0, 1.0);
            
            term_y = cos(2*pi * (y_prime / fibre_sep_y + phasefield_strength * (phasefield1 - 0.5)));
            term_z = cos(2*pi * (z_prime / fibre_sep_z + phasefield_strength * (phasefield2 - 0.5)));
            
            F = 0.5 + 0.5 * term_y .* term_z;
            F = F.^15;
            
            combined_noise = (1 - fibreness + fibreness * F) .* O_b;

        else
            % =============================================================
            % GLOBAL FALLBACK (For idealized 3D meshes without local fiber data)
            % =============================================================
            R = [ cos(phi) * cos(theta), sin(phi), -cos(phi) * sin(theta);
                 -cos(theta) * sin(phi), cos(phi),  sin(phi) * sin(theta);
                  sin(theta)           , 0       ,  cos(theta) ];

            P_rot = R * points'; 
            
            P_f = [ P_rot(1,:) / align_y^(1/3) / align_z^(1/3);
                    P_rot(2,:) * align_y^(2/3) / align_z^(1/3);
                    P_rot(3,:) / align_y^(1/3) * align_z^(2/3) ] / feature_size;

            O_b = Octave3D(P_f, 4, roughness, perm_table, offsets1);

            n_fibres_similarity = 4; wiggle_feature_length = 4; phasefield_strength = 5;
            phase_pts = [ P_rot(1,:) / wiggle_feature_length;
                          P_rot(2,:) / (n_fibres_similarity * fibre_sep_y);
                          P_rot(3,:) / (n_fibres_similarity * fibre_sep_z) ];
            
            phasefield1 = Octave3D(phase_pts, 4, 0.5, perm_table2, offsets1);
            phasefield2 = Octave3D(phase_pts, 4, 0.5, perm_table3, offsets2); 
            
            term_y = cos(2*pi * (P_rot(2,:) / fibre_sep_y + phasefield_strength * (phasefield1 - 0.5)));
            term_z = cos(2*pi * (P_rot(3,:) / fibre_sep_z + phasefield_strength * (phasefield2 - 0.5)));
            
            F = 0.5 + 0.5 * term_y .* term_z;
            F = F.^15;
            
            combined_noise = (1 - fibreness + fibreness * F) .* O_b;
        end

        % STEP 3: DENSITY VARIATION FIELD (O_d)
        O_d = Octave3D(points' / patch_size, 4, 0.5, perm_table4, offsets1);
        
        combined_noise = combined_noise + patchiness * O_d;

        noise_comps.O_b = O_b;
        noise_comps.F = F;
        noise_comps.O_d = O_d;
        noise_comps.G = combined_noise;

    % =====================================================================
    % BRANCH: 2D GENERATION
    % =====================================================================
    else
        % Extract 2D Params
        fibre_sep   = params.fibre_sep;
        align       = params.alignment;
        orientation = params.orientation;

        % Handle Non-Fibrous Case
        if isnan(fibreness)
            fibreness = 0; fibre_sep = 1;
        end

        points = mesh.mm_points(:, 1:2); % Always N x 2
        
        P_base = points' / feature_size; % Transpose to 2xN format for MEX engine
        
        if isfield(mesh, 'fibers')
            % Extract and normalize patient fibers
            f0 = mesh.fibers(:, 1:2);
            f0 = f0 ./ max(sqrt(sum(f0.^2, 2)), 1e-8);
            
            % Apply global orientation angle as a rotational offset
            R_off = [cos(orientation), -sin(orientation); 
                     sin(orientation),  cos(orientation)];
            f0_rot = (R_off * f0')'; 
            
            % STEP 1: MAIN FIBROSIS NOISE FIELD (O_b)
            O_b = Octave2D(P_base, 4, roughness, perm_table, offsets1, f0_rot', align);

            % STEP 2: FIBRE-SELECTING FIELD (F)
            centroid = mean(points, 1);
            v = points - centroid;
            t0 = [-f0_rot(:,2), f0_rot(:,1)]; 
            y_prime = sum(v .* t0, 2)'; 
            
            n_fibres_similarity = 4; wiggle_feature_length = 4; phasefield_strength = 5;
            P_veins = points' / wiggle_feature_length;
            
            phasefield = Octave2D(P_veins, 4, 0.5, perm_table2, offsets1, f0_rot', 1.0);
            
            F = 0.5 + 0.5 * cos(2*pi * (y_prime / fibre_sep + phasefield_strength * (phasefield - 0.5)));
            F = F.^15;
            
            combined_noise = (1 - fibreness + fibreness * F) .* O_b;
            
        else
            % =============================================================
            % GLOBAL FALLBACK (For idealized 2D meshes without local fiber data)
            % =============================================================
            
            % Simple global rotation based on user-provided orientation
            R = [ cos(orientation), sin(orientation);
                 -sin(orientation), cos(orientation) ]; 
            R_points = R * points';
            
            P_f = [ R_points(1,:) / sqrt(align);
                    R_points(2,:) * sqrt(align) ] / feature_size;
                    
            O_b = Octave2D(P_f, 4, roughness, perm_table, offsets1);
            
            n_fibres_similarity = 4; wiggle_feature_length = 4; phasefield_strength = 5;
            phase_pts = [ R_points(1,:) / wiggle_feature_length;
                          R_points(2,:) / (n_fibres_similarity * fibre_sep) ];
            phasefield = Octave2D(phase_pts, 4, 0.5, perm_table2, offsets1);
            
            F = 0.5 + 0.5 * cos(2*pi * (R_points(2,:) / fibre_sep + phasefield_strength * (phasefield - 0.5)));
            F = F.^15;
            combined_noise = (1 - fibreness + fibreness * F) .* O_b;
        end

        % STEP 3: DENSITY VARIATION FIELD (O_d)
        O_d = Octave2D(points' / patch_size, 4, 0.5, perm_table3, offsets1);
        
        combined_noise = combined_noise + patchiness * O_d;

        noise_comps.O_b = O_b;
        noise_comps.F = F;
        noise_comps.O_d = O_d;
        noise_comps.G = combined_noise;
    end
end