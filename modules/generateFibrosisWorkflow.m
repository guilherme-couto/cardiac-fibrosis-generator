function [output_data] = generateFibrosisWorkflow(domain_cfg, geom_cfg, bz_cfg, noise_params, run_cfg)
% GENERATEFIBROSISWORKFLOW Orchestrates the Geometry-Noise-Threshold (GNT) pipeline.
%
% This master function handles mesh construction/loading, geometrical masking, 
% Perlin noise evaluation, and morphological composition to generate synthetic 
% histology. It supports both analytical ideal grids and custom patient meshes.
%
% INPUTS:
%   domain_cfg   - Struct: Spatial dimensions (.Lx, .Ly, .Lz, .dx) or mesh file path.
%   geom_cfg     - Struct: Defines analytical boundaries (shape, center, width, depth).
%   bz_cfg       - Struct: Defines Border Zone (Gray Zone) topology and decay metric.
%   noise_params - Struct: Procedural Perlin variables (fibreness, patchiness, etc.).
%   run_cfg      - Struct: Execution modes ('single-pass'/'composition') and output flags.
%
% OUTPUTS:
%   output_data  - Struct containing the generated data:
%       .presence  - [N x 1] or [Ny x Nx x Nz] Boolean array of collagen placement.
%       .bz_layers - [N x 1] or [Ny x Nx x Nz] Integer array of topological layers.
%       .mesh      - Struct containing the processed spatial point-cloud.
%       .params    - Struct of the original noise parameters used.
%       .is3D      - Boolean flag mapping dimensionality.

    t_workflow = tic;
    fprintf('\n====  FIBROSIS GENERATION FRAMEWORK (%s - %s mode)  ====\n\n', run_cfg.fibrosis_type, run_cfg.mode);

    %% 1. BUILD OR LOAD MESH
    t_build_load_mesh = tic;
    if strcmp(domain_cfg.dimension_mode, 'CUSTOM')
        fprintf('[+] Loading patient-specific mesh from: %s\n', domain_cfg.mesh_file);
        mesh = readCustomMesh(domain_cfg.mesh_file); 

        min_b = min(mesh.points, [], 1);
        max_b = max(mesh.points, [], 1);
        max_range = max(max_b - min_b);
        
        % Forcing cm scale (The original generator was calibrated for mm)
        if max_range > 1000
            scale_to_mm = 1e-4; 
            fprintf('      |-> [!] Micrometer coordinate system detected. Scaling to cm.\n');
        elseif max_range < 10
            scale_to_mm = 1e2;  
            fprintf('      |-> [!] Meter coordinate system detected. Scaling to cm.\n');
        else
            scale_to_mm = 1e-1; 
            fprintf('      |-> [!] Millimeter coordinate system detected. Scaling to cm.\n');
        end
        
        mesh.points = mesh.points * scale_to_mm;
        mesh.Lx = (max_b(1) - min_b(1)) * scale_to_mm;
        mesh.Ly = (max_b(2) - min_b(2)) * scale_to_mm;
        mesh.Lz = (max_b(3) - min_b(3)) * scale_to_mm;
        if mesh.Lz == 0, mesh.Lz = 1e-6; end

        fprintf('      |-> [✓] Mesh ready with %d points. Dimensions: [%.2f x %.2f x %.2f] cm.\n', ...
            mesh.num_points, mesh.Lx, mesh.Ly, mesh.Lz);
    else
        fprintf('[+] Building analytical structured grid\n');
        mesh = buildMesh(domain_cfg);
    end
    elapsed_build_load_mesh = toc(t_build_load_mesh);

    %% 2. & 3. COMPUTE GEOMETRY & GENERATE FIBROSIS PATTERN
    t_pattern_gen = tic;
    if strcmp(domain_cfg.dimension_mode, 'CUSTOM')
        % CUSTOM MESH LOGIC
        fprintf('[+] Defining custom geometry\n');
        
        core_mask = (mesh.tags == 1 | mesh.tags == 2);
        core_indices_all = find(core_mask);
        
        if isempty(core_indices_all)
            error('      |-> [x] Generation halted: No fibrotic regions (Tags 1 or 2) found in custom mesh.');
        end
        
        sub_mesh = mesh;
        sub_mesh.points = mesh.points(core_indices_all, :);
        sub_mesh.num_points = length(core_indices_all);
        sub_mesh.tags = mesh.tags(core_indices_all);
        
        if isfield(mesh, 'fibers'), sub_mesh.fibers = mesh.fibers(core_indices_all, :); end
        if isfield(mesh, 'sheet'),  sub_mesh.sheet  = mesh.sheet(core_indices_all, :);  end
        if isfield(mesh, 'normal'), sub_mesh.normal = mesh.normal(core_indices_all, :); end

        geometry.core_mask = true(sub_mesh.num_points, 1);
        geometry.core_indices = (1:sub_mesh.num_points)';
        geometry.num_core_points = sub_mesh.num_points;
        geometry.is3D = mesh.is3D;
        geometry.tags = sub_mesh.tags;
        geometry.bz_map = zeros(sub_mesh.num_points, 1);
        geometry.total_layers = 0;
        
        if ~isempty(mesh.gz_weights)
            geometry.weights = mesh.gz_weights(core_indices_all);
        else
            geometry.weights = []; 
        end
        
        geometry.true_core_mask = (sub_mesh.tags == 1);
        geometry.num_true_core = sum(geometry.true_core_mask);
        if geometry.num_true_core == 0
            geometry.true_core_mask = (sub_mesh.tags == 2);
            geometry.num_true_core = sum(geometry.true_core_mask);
        end
        geometry.true_core_indices = find(geometry.true_core_mask);

        %% 3. GENERATE PATTERN (CUSTOM)
        fprintf('[+] Generating fibrosis pattern\n');
        
        if strcmpi(run_cfg.mode, 'single-pass')
            [combined_noise, ~] = computeNoiseField(sub_mesh, noise_params);
            local_res = applyFibrosisThreshold(combined_noise, geometry, run_cfg.target_density, bz_cfg);
        elseif strcmpi(run_cfg.mode, 'composition')
            local_res = runCompositionLoop(sub_mesh, geometry, noise_params, run_cfg.target_density, bz_cfg);
        else
            error('      |-> [x] Unknown generation mode: %s', run_cfg.mode);
        end
        
        final_result.presence = false(mesh.num_points, 1);
        final_result.presence(core_indices_all) = local_res.presence;
        
        final_result.bz_layers = -ones(mesh.num_points, 1);
        final_result.bz_layers(core_indices_all) = local_res.bz_layers;
        final_result.fc_density = local_res.fc_density; 
        
        core_presence = final_result.presence(mesh.tags == 1);
        gz_presence = final_result.presence(mesh.tags == 2);
        
        num_core_total = length(core_presence);
        num_core_col = sum(core_presence);
        core_pct = (num_core_col / max(1, num_core_total)) * 100;
        
        num_gz_total = length(gz_presence);
        num_gz_col = sum(gz_presence);
        gz_pct = (num_gz_col / max(1, num_gz_total)) * 100;

        fprintf('      |-> [✓] Synthesis complete. Collagen distribution:\n');
        fprintf('               |-> Fibrotic Core: %d / %d (%.2f%%)\n', num_core_col, num_core_total, core_pct);
        fprintf('               |-> Gray Zone:     %d / %d (%.2f%%)\n', num_gz_col, num_gz_total, gz_pct);

    else
        % IDEALIZED GEOMETRY LOGIC
        fprintf('[+] Defining idealized geometry\n');
        geometry = computeFibrosisGeometry(mesh, geom_cfg, bz_cfg);
        
        geometry.true_core_mask = true(geometry.num_core_points, 1);
        geometry.num_true_core = geometry.num_core_points;
        geometry.true_core_indices = geometry.core_indices;

        if geometry.num_core_points == 0
            error('      |-> [x] Generation halted: Idealized core dimensions resulted in an empty region.');
        end
        
        %% 3. GENERATE PATTERN (Idealized)
        fprintf('[+] Generating fibrosis pattern\n');
        if strcmpi(run_cfg.mode, 'uniform')
            final_result = generateUniformPattern(geometry, run_cfg.target_density, noise_params.seed, bz_cfg);
        elseif strcmpi(run_cfg.mode, 'single-pass')
            [combined_noise, ~] = computeNoiseField(mesh, noise_params);
            final_result = applyFibrosisThreshold(combined_noise, geometry, run_cfg.target_density, bz_cfg);
        elseif strcmpi(run_cfg.mode, 'composition')
            final_result = runCompositionLoop(mesh, geometry, noise_params, run_cfg.target_density, bz_cfg);
        end
    end
    elapsed_pattern_gen = toc(t_pattern_gen);

    %% 4. PACK AND EXPORT RESULTS
    output_data.presence = final_result.presence;
    output_data.bz_layers = final_result.bz_layers;
    output_data.mesh = mesh;
    output_data.params = noise_params;
    output_data.is3D = mesh.is3D;
    
    if strcmpi(domain_cfg.dimension_mode, 'CUSTOM')
        if run_cfg.save_mesh
            t_export = tic;
            fprintf('[+] Exporting computational meshes\n');
            writeCustomMesh(mesh, final_result.presence, run_cfg.output_name);
            elapsed_export = toc(t_export);
        end
        if run_cfg.save_fig & !output_data.is3D
            t_visualization = tic;
            fprintf('[+] Generating visualizations\n');
            saveCustomPatternVisualisation(mesh, run_cfg.output_name, run_cfg);
            elapsed_visualization = toc(t_visualization);
        end
    else
        % Analytical Mode
        if run_cfg.save_mesh
            t_export = tic;
            fprintf('[+] Exporting computational meshes\n');
            writeALGFromPatternData(output_data, run_cfg.output_name);
            elapsed_export = toc(t_export);
        end
        if run_cfg.save_fig
            t_visualization = tic;
            fprintf('[+] Generating visualizations\n');
            savePatternVisualisation(output_data, final_result, run_cfg);
            elapsed_visualization = toc(t_visualization);
        end
    end

    elapsed_workflow = toc(t_workflow);

    fprintf('[✓] Total workflow completed in %.2f seconds.\n', elapsed_workflow);
    fprintf('      |-> Mesh construction/loading completed in %.2f seconds.\n', elapsed_build_load_mesh);
    fprintf('      |-> Pattern generation completed in %.2f seconds.\n', elapsed_pattern_gen);
    if run_cfg.save_mesh
        fprintf('      |-> Export completed in %.2f seconds.\n', elapsed_export);
    end
    if run_cfg.save_fig
        fprintf('      |-> Visualization generated in %.2f seconds.\n', elapsed_visualization);
    end

    fprintf('\n====  FRAMEWORK COMPLETE  ====\n\n');
end