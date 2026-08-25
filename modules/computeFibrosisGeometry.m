function geometry = computeFibrosisGeometry(mesh, geom_cfg, bz_cfg)
% COMPUTEFIBROSISGEOMETRY Defines spatial boundaries for fibrosis generation.
%
% Creates logical masking arrays that dictate where the fibrosis core and 
% associated border zone layers are located. Supports analytical 2D/3D shapes 
% embedded in structured grids, and automatically maps unstructured patient tags.
%
% INPUTS:
%   mesh     - Struct containing point cloud data (.points, .tags, dimensions).
%   geom_cfg - Struct defining shape properties (shape type, center, width, depth).
%   bz_cfg   - Struct defining Border Zone rules (distance metric, max layers).
%
% OUTPUTS:
%   geometry - Struct containing:
%       .core_mask       - Logical array indicating valid core target areas.
%       .core_indices    - Numerical indices of the core_mask.
%       .num_core_points - Total number of valid elements.
%       .bz_map          - Array classifying elements into discrete spatial layers.
%       .total_layers    - Maximum layer depth computed by the distance transform.

    % BRANCH 1: CUSTOM MESH (PATIENT GEOMETRY)
    if isfield(mesh, 'tags')
        % In custom mode, geometry is dictated by the patient's mesh tags, not by analytical shapes.
        % User rule: Tags 1 and 2 are considered Fibrotic Core. 0 is Healthy.
        core_mask = (mesh.tags == 1 | mesh.tags == 2);
        
        geometry.core_mask = core_mask;       % N x 1 Logical Vector
        geometry.core_indices = find(core_mask);
        geometry.num_core_points = length(geometry.core_indices);
        
        % Border Zone logic for custom meshes is associated with native tags (gray zone) and not computed via distance transforms
        final_bz_map = -ones(mesh.num_points, 1);
        final_bz_map(core_mask) = 0;
        
        geometry.bz_map = final_bz_map;
        geometry.total_layers = 0;
        geometry.is3D = mesh.is3D;
        
        fprintf('  -> Custom Geometry: Identified %d fibrotic core elements.\n', geometry.num_core_points);
        return;
    end

    % BRANCH 2: ANALYTICAL GEOMETRY (IDEAL 2D/3D GRIDS)
    Nx = mesh.Nx;
    Ny = mesh.Ny;
    Nz = mesh.Nz;
    is3D = mesh.is3D;
    
    % Initialize Core Mask with appropriate dimensions
    if is3D
        core_mask = false(Ny, Nx, Nz);
    else
        core_mask = false(Ny, Nx);
    end

    % Dimensions of the fibrotic core in PIXELS for embedded shapes
    if isfield(geom_cfg, 'fc_width'), w_px = geom_cfg.fc_width / mesh.dx; end
    if isfield(geom_cfg, 'fc_height'), h_px = geom_cfg.fc_height / mesh.dx; end
    if isfield(geom_cfg, 'fc_center_x'), cx_px = geom_cfg.fc_center_x / mesh.dx; end
    if isfield(geom_cfg, 'fc_center_y'), cy_px = geom_cfg.fc_center_y / mesh.dx; end
    
    % DEFINE CORE SHAPE
    switch lower(geom_cfg.shape)
        
        case 'full'
            % Works for both 2D and 3D
            if is3D
                core_mask = true(Ny, Nx, Nz);
            else
                core_mask = true(Ny, Nx);
            end

        % 2D SHAPES
        case 'rectangle'
            if is3D, error('Shape "rectangle" is for 2D meshes only. Use "box" for 3D.'); end
            
            [X, Y] = meshgrid(1:Nx, 1:Ny);
            x_min = cx_px - w_px/2; x_max = cx_px + w_px/2;
            y_min = cy_px - h_px/2; y_max = cy_px + h_px/2;
            
            core_mask = (X >= x_min & X <= x_max & Y >= y_min & Y <= y_max);
            
        case 'ellipse'
            if is3D, error('Shape "ellipse" is for 2D meshes only. Use "ellipsoid" for 3D.'); end
            
            [X, Y] = meshgrid(1:Nx, 1:Ny);
            a = w_px/2; 
            b = h_px/2;
            
            % Normalized ellipse equation: ((x-cx)/a)^2 + ((y-cy)/b)^2 <= 1
            dist_sq = ((X - cx_px).^2) / a^2 + ((Y - cy_px).^2) / b^2;
            core_mask = (dist_sq <= 1);

        % 3D SHAPES
        case 'box'
            if ~is3D, error('Shape "box" requires a 3D mesh (Lz > 0).'); end
            
            % Z-Dimension handling
            d_px = getDepthPx(geom_cfg, mesh, Nz);
            cz_px = getCenterZPx(geom_cfg, mesh, Nz);

            [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
            
            x_min = cx_px - w_px/2; x_max = cx_px + w_px/2;
            y_min = cy_px - h_px/2; y_max = cy_px + h_px/2;
            z_min = cz_px - d_px/2; z_max = cz_px + d_px/2;
            
            core_mask = (X >= x_min & X <= x_max & ...
                         Y >= y_min & Y <= y_max & ...
                         Z >= z_min & Z <= z_max);

        case 'ellipsoid'
            if ~is3D, error('Shape "ellipsoid" requires a 3D mesh (Lz > 0).'); end
            
            % Z-Dimension handling
            d_px = getDepthPx(geom_cfg, mesh, Nz);
            cz_px = getCenterZPx(geom_cfg, mesh, Nz);

            [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
            
            a = w_px/2; 
            b = h_px/2;
            c = d_px/2;
            if c == 0, c = 1e-6; end % Prevent division by zero
            
            % Normalized ellipsoid equation: ((x-cx)/a)^2 + ((y-cy)/b)^2 + ((z-cz)/c)^2 <= 1
            dist_sq = ((X - cx_px).^2) / a^2 + ...
                      ((Y - cy_px).^2) / b^2 + ...
                      ((Z - cz_px).^2) / c^2;
            core_mask = (dist_sq <= 1);
            
        case 'stain'
            error('Stain shape logic needs migration.');
            
        otherwise
            error('Unknown shape: %s', geom_cfg.shape);
    end
    
    % BORDER ZONE CALCULATION
    % Using Distance Transform
    if bz_cfg.active
        % bwdist calculates distance to nearest non-zero pixel.
        % We want distance from the core (where core is 1).
        % So we compute distance on the INVERSE of the core.
        % 'chessboard' = 8-connectivity.
        % 'euclidean' = smooth circles.
        dist_map = bwdist(core_mask, bz_cfg.metric); 
        
        % Convert distance to integer layers (discretize continuous distances)
        % Distance 0 = Inside Core
        % Distance 1 = Immediate neighbor (Layer 1)
        layer_map = ceil(double(dist_map));
        
        % Cap at max layers
        if isfield(bz_cfg, 'max_layers') && ~isinf(bz_cfg.max_layers)
            mask_valid_bz = (layer_map > 0) & (layer_map <= bz_cfg.max_layers);
        else
            mask_valid_bz = (layer_map > 0);
        end
        
        % Create Final Maps
        % Convention: 0 = Core, 1..N = BZ, -1 = Healthy
        final_bz_map = -ones(size(core_mask));
        final_bz_map(core_mask) = 0;
        final_bz_map(mask_valid_bz) = layer_map(mask_valid_bz);
        
        total_layers = max(final_bz_map(:));
    else
        % No BZ logic
        final_bz_map = -ones(size(core_mask));
        final_bz_map(core_mask) = 0;
        total_layers = 0;
    end
    
    geometry.core_mask = core_mask;       % Logical Matrix
    geometry.core_indices = find(core_mask);
    geometry.num_core_points = length(geometry.core_indices);
    geometry.bz_map = final_bz_map;       % Matrix (-1, 0, 1..N)
    geometry.total_layers = total_layers;
    geometry.is3D = is3D;
end

% Helper Functions for 3D Params
function d_px = getDepthPx(geom_cfg, mesh, Nz)
    if isfield(geom_cfg, 'fc_depth')
        d_px = geom_cfg.fc_depth / mesh.dx;
    else
        d_px = Nz; % Default to full depth
    end
end

function cz_px = getCenterZPx(geom_cfg, mesh, Nz)
    if isfield(geom_cfg, 'fc_center_z')
        cz_px = geom_cfg.fc_center_z / mesh.dx;
    else
        cz_px = Nz / 2; % Default to center
    end
end
