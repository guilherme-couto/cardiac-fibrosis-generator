function mesh = buildMesh(domain_config)
% BUILDMESH Constructs an analytical Cartesian grid for procedural generation.
%
% Generates a structured point-cloud based on physical dimensions. Supports
% both 2D and 3D generation. The points are centered within each computational pixel/voxel.
%
% INPUTS:
%   domain_config - A struct containing:
%       .Lx : Width of the domain (cm)
%       .Ly : Height of the domain (cm)
%       .Lz : Depth of the domain (cm) [Optional. Defaults to 0 for 2D]
%       .dx : Spatial resolution (cm)
%
% OUTPUTS:
%   mesh - A standardized framework struct containing:
%       .points : [N x 3] Array of spatial coordinates
%       .Nx, .Ny, .Nz : Number of elements in each dimension
%       .dx : Spatial resolution
%       .is3D : Boolean flag indicating dimensionality

    Lx = domain_config.Lx;
    Ly = domain_config.Ly;
    dx = domain_config.dx;
    Lz = 0;
    
    if isfield(domain_config, 'Lz')
        Lz = domain_config.Lz;
    end

    Nx = ceil(Lx/dx);
    Ny = ceil(Ly/dx);

    if Lz == 0
        % 2D Generation
        xv = linspace(dx/2, dx*(Nx - 0.5), Nx);
        yv = linspace(dx/2, dx*(Ny - 0.5), Ny);
        [X,Y] = meshgrid(xv,yv);
        points = [X(:), Y(:)];
        mesh.Nz = 1;
        mesh.is3D = false;
    else
        % 3D Generation
        Nz = ceil(Lz/dx);
        xv = linspace(dx/2, dx*(Nx - 0.5), Nx);
        yv = linspace(dx/2, dx*(Ny - 0.5), Ny);
        zv = linspace(dx/2, dx*(Nz - 0.5), Nz);
        [X,Y,Z] = meshgrid(xv,yv,zv);
        points = [X(:), Y(:), Z(:)];
        mesh.Nz = Nz;
        mesh.is3D = true;
    end

    mesh.points = points;
    mesh.Nx = Nx;
    mesh.Ny = Ny;
    mesh.dx = dx;

    fprintf('      |-> Mesh ready with %d points. Dimensions: [%.2f x %.2f x %.2f] cm.\n', ...
            size(points, 1), Lx, Ly, Lz);
end