function mesh = readCustomMesh(filename)
% READCUSTOMMESH Parses and loads large-scale unstructured patient meshes.
%
% INPUTS:
%   filename - String containing the relative or absolute path to the mesh file.
%
% OUTPUTS:
%   mesh     - Framework-standardized struct containing:
%       .points        - [N x 3] Array of spatial coordinates.
%       .tags          - [N x 1] Array of anatomical tissue tags.
%       .num_points    - Total number of elements.
%       .is3D          - Boolean indicating dimensionality.
%       .fibers        - [N x 3] (or 2) Biological fiber vectors.
%       .sheet         - [N x 3] Transverse sheet vectors (3D only).
%       .normal        - [N x 3] Sheet normal vectors (3D only).
%       .gz_weights    - [N x 1] Optional Laplacian weights for Gray Zone.
%       .type          - Extension identifier ('alg' or 'vtu').
%       .original_file - Path to the source file for later injection.

    % =====================================================================
    % DEVELOPER CONFIGURATIONS
    % =====================================================================
    % Toggle to save/load parsed meshes as binary .mat files.
    % Saves minutes on consecutive runs with the same 3D grid.
    ENABLE_CACHE = true;
    CACHE_DIR = 'cached_meshes';

    [~, mesh_name, ext] = fileparts(filename);
    fprintf('  [i] Loading custom patient mesh: %s\n', filename);

    % =====================================================================
    % CACHE SYSTEM (Ultra-fast I/O bypass)
    % =====================================================================
    if ENABLE_CACHE
        if ~exist(CACHE_DIR, 'dir')
            mkdir(CACHE_DIR);
        end
        cache_file = fullfile(CACHE_DIR, sprintf('%s_cache.mat', mesh_name));

        if exist(cache_file, 'file')
            fprintf('  [i] Found pre-compiled mesh data in cached_meshes: %s. Loading...\n', cache_file);
            tic();
            loaded_data = load(cache_file, 'mesh');
            mesh = loaded_data.mesh;
            fprintf('  [✓] Mesh loaded from cache (took %.2f seconds).\n', toc());
            return;
        end
    end

    if strcmpi(ext, '.alg')
        fid = fopen(filename, 'r');
        if fid == -1
            error('  [x] Could not open file %s', filename);
        end

        tic();
        % 1. Peek at the first line to dynamically count columns
        tline = fgetl(fid);
        num_cols = length(strfind(tline, ',')) + 1;
        frewind(fid);
        elapsed_time = toc();
        fprintf('  [i] Detected %d columns in ALG file (took %.2f seconds)\n', num_cols, elapsed_time);

        if num_cols == 0
            fclose(fid);
            error('  [x] ALG file appears to be empty or corrupted.');
        end

        % 2. Smart Format Specifier (Skip unused columns to save RAM and CPU)
        tic();
        fmt = '';
        col_map = zeros(num_cols, 1);
        current_idx = 1;

        for c = 1:num_cols
            % Define which columns actually matter to the framework
            is_useful = (c <= 3) || (c == 7) || (c >= 9 && c <= 18);
            
            if is_useful
                fmt = [fmt, '%f '];
                col_map(c) = current_idx;
                current_idx = current_idx + 1;
            else
                % %*f tells textscan to read but completely discard the data
                fmt = [fmt, '%*f '];
            end
        end
        elapsed_time = toc();
        fprintf('  [i] Prepared selective format specifier (took %.2f seconds)\n', elapsed_time);

        % 3. Native C++ Block Read (No cell2mat duplication)
        tic();
        % Removed CommentStyle to unlock maximum parsing speed.
        % CollectOutput automatically concatenates the results into a single matrix.
        data_cell = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
        fclose(fid);
        
        alg_data = data_cell{1};
        elapsed_time = toc();
        fprintf('  [i] Parsed and allocated directly to matrix (took %.2f seconds)\n', elapsed_time);
        
        % =================================================================
        % MESH EXTRACTIONS (Using the dynamic column map)
        % =================================================================
        mesh.points = alg_data(:, col_map(1:3));
        mesh.tags = alg_data(:, col_map(7));
        mesh.num_points = size(alg_data, 1);
        mesh.is3D = (num_cols >= 18);
        
        % Gray Zone Extraction
        if num_cols >= 18
            mesh.gz_weights = alg_data(:, col_map(18));
        elseif num_cols == 12
            mesh.gz_weights = alg_data(:, col_map(12));
        else
            mesh.gz_weights = [];
        end

        mesh.type = 'alg';
        mesh.original_file = filename; 

        % Biological Tensors Extraction
        if num_cols >= 11 && num_cols < 18
            mesh.fibers = alg_data(:, col_map(9:11));
        elseif num_cols >= 18
            mesh.fibers = alg_data(:, col_map(9:11));
            mesh.sheet  = alg_data(:, col_map(12:14));
            mesh.normal = alg_data(:, col_map(15:17));
        end

    elseif strcmpi(ext, '.vtu')
        % Delegate binary XML parsing and centroid calculation to Python
        temp_mat = 'temp_vtu_in.mat';
        
        bridge_script = fullfile(pwd, 'helpers', 'vtu_bridge.py');
        cmd = sprintf("python3 '%s' read '%s' '%s'", bridge_script, filename, temp_mat);
        
        [status, cmd_out] = system(cmd);
        if status ~= 0
            error('  [x] Failed to read VTU via Python:\n%s', cmd_out);
        end
        
        % Load parsed data
        load(temp_mat, 'points', 'tags');
        mesh.points = points;
        mesh.tags = double(tags);
        
        mesh.type = 'vtu';
        mesh.original_file = filename; 
        
        % Clean up
        delete(temp_mat);
    else
        error('  [x] Unsupported mesh extension: %s. Use .alg or .vtu.', ext);
    end
    
    % Core standardizations for the framework
    mesh.num_points = length(mesh.tags);
    
    % Check if 3D (If all Z coordinates are the same, it's a 2D slice)
    z_coords = mesh.points(:, 3);
    mesh.is3D = (max(z_coords) - min(z_coords)) > 1e-6;
    
    mesh.min_bounds = min(mesh.points, [], 1);
    mesh.max_bounds = max(mesh.points, [], 1);
    
    fprintf('  [✓] Loaded %d elements. 3D Mode: %d\n', mesh.num_points, mesh.is3D);
    fprintf('       |-> Tag 0 (Healthy): %d elements\n', sum(mesh.tags == 0));
    fprintf('       |-> Tag 1 (Core):    %d elements\n', sum(mesh.tags == 1));
    fprintf('       |-> Tag 2 (GZ):      %d elements\n', sum(mesh.tags == 2));

    % =====================================================================
    % SAVE CACHE (For future runs)
    % =====================================================================
    if ENABLE_CACHE
        fprintf('  [i] Saving uncompressed parsed mesh to %s...\n', cache_file);
        save('-binary', cache_file, 'mesh');
    end
end