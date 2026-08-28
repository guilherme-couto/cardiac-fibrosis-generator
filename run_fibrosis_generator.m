function [output_data] = run_fibrosis_generator(fibrosis_type, density, seed, angle_deg_or_vec, dimension_mode, domain_dims_or_path, shape, core_dims, output_filename, save_mesh, save_figure)
% RUN_FIBROSIS_GENERATOR Entry point for the fibrosis generation framework.
%
% This function acts as a wrapper to sanitize user inputs from scripts or 
% CLI arguments before calling the core generation workflow. 
%
% PREREQUISITES:
%   The 'image' package must be installed in Octave (pkg install -forge image).
%
% INPUTS:
%   fibrosis_type       : String ('compact', 'interstitial', 'diffuse', 'patchy', 'uniform')
%   density             : Float [0, 1] Target volumetric density.
%   seed                : Integer RNG seed for exact reproducibility.
%   angle_deg_or_vec    : Float (2D orientation) or Vector [Phi, Theta] for 3D fiber angles. Not needed for CUSTOM meshes.
%   dimension_mode      : String ('2D', '3D', 'CUSTOM').
%   domain_dims_or_path : Array [dx, Lx, Ly, Lz] in cm. For CUSTOM, pass the mesh file path string.
%   shape               : String ('full', 'ellipse', 'rectangle', 'box', 'ellipsoid'). Not needed for CUSTOM meshes.
%   core_dims           : Array [width, height, depth] in cm. Ignored if shape is 'full'.
%   output_filename     : String base name for the output files (no extensions).
%   save_mesh           : Boolean flag to write the physical arrays to disk.
%   save_figure         : Boolean flag to render PNG/GIF visualizations.

    %% 0. FRAMEWORK SETUP
    addpath(pwd);
    addpath(fullfile(pwd, 'modules'));
    addpath(fullfile(pwd, 'helpers'));

    try 
        pkg load image; 
    catch 
        error('[x] Octave "image" package not found. Please run: pkg install -forge image'); 
    end
    more off; 

    %% 1. INPUT SANITIZATION
    if isempty(fibrosis_type) || ~ischar(fibrosis_type), error('[x] fibrosis_type must be a non-empty string.'); end
    if isempty(density) || ~isnumeric(density) || density < 0 || density > 1, error('[x] density must be a numeric value in [0, 1].'); end
    if isempty(seed) || ~isnumeric(seed) || seed < 0, error('[x] seed must be a non-negative integer.'); end
    if isempty(angle_deg_or_vec) && ~strcmpi(upper(dimension_mode), 'CUSTOM'), error('[x] angle_deg_or_vec must be provided for 2D/3D meshes.'); end
    if isempty(dimension_mode) || ~ischar(dimension_mode), error('[x] dimension_mode must be a non-empty string.'); end
    if isempty(domain_dims_or_path), error('[x] domain_dims_or_path must be provided.'); end
    if (isempty(shape) || ~ischar(shape)) && ~strcmpi(upper(dimension_mode), 'CUSTOM'), error('[x] shape must be provided for 2D/3D meshes and should be a string.'); end
    if isempty(core_dims) && ~strcmpi(shape, 'full') && ~strcmpi(upper(dimension_mode), 'CUSTOM'), error('[x] core_dims must be provided for non-full shapes.'); end
    if isempty(output_filename) || ~ischar(output_filename), error('[x] output_filename must be a non-empty string.'); end
    if isempty(save_mesh) || ~islogical(save_mesh), error('[x] save_mesh must be a boolean flag.'); end
    if isempty(save_figure) || ~islogical(save_figure), error('[x] save_figure must be a boolean flag.'); end

    dimension_mode = upper(char(dimension_mode)); 

    %% 2. FIBROSIS DATABASE DEFINITION
    override_mode = '';
    raw_type = lower(fibrosis_type);
    
    % If the user appends '_single' or '_comp' to the fibrosis type, override the default mode
    if endsWith(raw_type, '_single')
        override_mode = 'single-pass';
        raw_type = strrep(raw_type, '_single', '');
    elseif endsWith(raw_type, '_comp')
        override_mode = 'composition';
        raw_type = strrep(raw_type, '_comp', '');
    end
    
    % Core Parameters: [fibreness, fibre_sep, patchiness, feature_size, roughness, patch_size, alignment, 0]
    % Obs: Values are from Lawson et al. original paper
    switch raw_type
        case 'compact',      base_params = [NaN, NaN, 0.44, 0.96, 0.59, 2.03, 2.47, 0]; mode = 'single-pass';
        case 'interstitial', base_params = [0.30, 0.31, 0.32, 0.24, 0.96, 4.67, 1.89, 0]; mode = 'composition';
        case 'diffuse',      base_params = [NaN, NaN, 0.49, 0.07, 0.44, 1.22, 2.17, 0]; mode = 'composition';
        case 'patchy',       base_params = [0.38, 0.31, 0.42, 0.32, 0.78, 2.10, 2.50, 0]; mode = 'composition';
        case 'uniform',      base_params = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0]; mode = 'uniform';
        otherwise, error('  [x] Unknown fibrosis phenotype: %s', raw_type);
    end
    
    if ~isempty(override_mode)
        mode = override_mode;
    end
    fibrosis_type = raw_type;

    %% 3. BUILD CONFIGURATION STRUCTS
    domain_config.dimension_mode = dimension_mode; 
    
    if strcmp(dimension_mode, 'CUSTOM')
        angle_deg_or_vec = [];
        core_dims = [];
        shape = [];
        if ~ischar(domain_dims_or_path)
            error('[x] For CUSTOM meshes, domain_dims_or_path must be a string path to the mesh file.');
        end
        if ~isfile(domain_dims_or_path)
            error('[x] The specified custom mesh file does not exist: %s', domain_dims_or_path);
        end
        domain_config.mesh_file = domain_dims_or_path;
        domain_config.dx = 1; domain_config.Lx = 1; domain_config.Ly = 1; domain_config.Lz = 1;
        is3D = false; 
    else
        angle_deg_or_vec = parseVectorInput(angle_deg_or_vec);
        core_dims = parseVectorInput(core_dims);
        domain_dims_or_path = parseVectorInput(domain_dims_or_path);
        if length(domain_dims_or_path) < 3
            error('[x] domain_dims_or_path requires at least 3 elements: [dx, Lx, Ly]');
        end
        domain_config.dx = domain_dims_or_path(1);
        domain_config.Lx = domain_dims_or_path(2);
        domain_config.Ly = domain_dims_or_path(3);
        
        if length(domain_dims_or_path) == 4
            domain_config.Lz = domain_dims_or_path(4);
        else
            domain_config.Lz = 0.0;
        end

        if strcmp(dimension_mode, '3D')
            if domain_config.Lz <= 0
                warning('  [WARNING] 3D mode selected but Lz=0. Forcing Lz=dx.'); 
                domain_config.Lz = domain_config.dx; 
            end
            is3D = true;
        elseif strcmp(dimension_mode, '2D')
            domain_config.Lz = 0.0; 
            is3D = false;
        end
    end

    geom_config.shape = shape;

    if ~strcmpi(shape, 'full') && ~strcmpi(dimension_mode, 'CUSTOM')
        if length(core_dims) < 2
            error('[x] core_dims must define at least [width, height]');
        end

        geom_config.fc_width  = core_dims(1); 
        geom_config.fc_height = core_dims(2); 
        if length(core_dims) == 3
            geom_config.fc_depth  = core_dims(3); 
        end

        if geom_config.fc_width > domain_config.Lx || geom_config.fc_height > domain_config.Ly
            error('[x] Fibrotic core dimensions exceed analytical domain limits.');
        end

        geom_config.fc_center_x = domain_config.Lx / 2;
        geom_config.fc_center_y = domain_config.Ly / 2;
        geom_config.fc_center_z = domain_config.Lz / 2;

        if is3D && any(strcmpi(shape, {'rectangle', 'ellipse'}))
            error('  [x] Shape "%s" is strictly 2D. Use "box" or "ellipsoid".', shape);
        elseif ~is3D && any(strcmpi(shape, {'box', 'ellipsoid'}))
            error('  [x] Shape "%s" is strictly 3D. Use "rectangle" or "ellipse".', shape);
        end
    end
    
    if strcmpi(shape, 'full') || strcmpi(dimension_mode, 'CUSTOM')
        bz_config.active = false;
        bz_config.max_layers = 0;
    else
        bz_config.active = true;

        proportional_bz_thickness = 0.1; % 10% of the smallest core dimension
        if is3D
            bz_thickness = proportional_bz_thickness * min([geom_config.fc_width, geom_config.fc_height, geom_config.fc_depth]);
        else
            bz_thickness = proportional_bz_thickness * min([geom_config.fc_width, geom_config.fc_height]);
        end
        bz_config.max_layers = ceil(bz_thickness / domain_config.dx);
    end

    noise_params = vectorToParams(base_params, angle_deg_or_vec, is3D);
    noise_params.seed = seed; 

    run_config.fibrosis_type = fibrosis_type;
    run_config.mode = mode;
    run_config.target_density = density;
    run_config.seed = seed;
    run_config.output_name = output_filename;
    run_config.save_mesh = save_mesh;
    run_config.save_fig = save_figure;

    [output_dir, ~, ~] = fileparts(output_filename);
    if ~isempty(output_dir)
        [~, ~, ~] = mkdir(output_dir);
    end

    try
        output_data = generateFibrosisWorkflow(domain_config, geom_config, bz_config, noise_params, run_config);
    catch ME
        error('\n  [CRITICAL FAILURE] %s\n  Traceback: %s (Line %d)\n', ME.message, ME.stack(1).name, ME.stack(1).line);
    end
end

% Helpers
function val = parseVectorInput(input_val)
    if ischar(input_val) || isstring(input_val)
        clean_str = strrep(strrep(input_val, '[', ''), ']', '');
        val = str2num(clean_str);
    else
        val = input_val;
    end
end

function params_struct = vectorToParams(vec, angle_deg_or_vec, is3D)
    if is3D && length(vec) == 8
        vec = convert_params_2D_to_3D(vec);
    end
    
    params_struct.fibreness     = vec(1);
    params_struct.fibre_sep_y   = vec(2);
    params_struct.patchiness    = vec(3);
    params_struct.feature_size  = vec(4);
    params_struct.roughness     = vec(5);
    params_struct.patch_size    = vec(6);
    params_struct.alignment_y   = vec(7);
    
    if is3D
        params_struct.fibre_sep_z   = vec(8);
        params_struct.alignment_z   = vec(9);
        if length(angle_deg_or_vec) == 2
            params_struct.phi           = deg2rad(angle_deg_or_vec(1));
            params_struct.theta         = deg2rad(angle_deg_or_vec(2));
        elseif isscalar(angle_deg_or_vec)
            params_struct.phi           = deg2rad(angle_deg_or_vec);
            params_struct.theta         = deg2rad(angle_deg_or_vec);
        else
            error('  [x] For 3D meshes, angle vector must be a scalar or a [Phi, Theta] pair.');
        end
        params_struct.orientation = 0; 
    else
        params_struct.fibre_sep_z   = 1;
        params_struct.alignment_z   = 1;
        params_struct.phi           = 0;
        params_struct.theta         = 0;
        params_struct.orientation = deg2rad(angle_deg_or_vec);
    end
    
    params_struct.fibre_sep = params_struct.fibre_sep_y;
    params_struct.alignment = params_struct.alignment_y;
end

function params3D = convert_params_2D_to_3D(params2D)
    fibreness     = params2D(1);
    fibre_sep_y   = params2D(2);
    patchiness    = params2D(3);
    feature_size  = params2D(4);
    roughness     = params2D(5);
    patch_size    = params2D(6);
    alignment_y   = params2D(7);

    K = 1; 

    fibre_sep_z = fibre_sep_y;
    feature_size = feature_size .* alignment_y.^(-1/6) * K^(-1/3);
    alignment_z = alignment_y * K;

    phi = 0;    
    theta = 0;  

    params3D = [fibreness, fibre_sep_y, patchiness, feature_size, roughness, patch_size, alignment_y, fibre_sep_z, alignment_z, phi, theta];
end