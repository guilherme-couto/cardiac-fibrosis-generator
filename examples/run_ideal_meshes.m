% EXAMPLE: ANALYTICAL GRIDS
% Demonstrates the creation of standard geometrical domains and morphological rendering (GIFs and PNGs).

clc;
fprintf("Cardiac Fibrosis Framework Example: Analytical Ideal Grids\n\n");

% --- SCENARIO 1: 2D Interstitial Ellipse ---
scn1.id = 3;
scn1.type = 'interstitial'; 
scn1.dim_mode = '2D'; 
scn1.shape = 'ellipse'; 
scn1.angle = 60;
scn1.density = 0.10;
scn1.domain = [0.01, 2.0, 4.0];  % [dx, Lx, Ly]
scn1.core   = [1.5, 3.0];        % [width, height]
scn1.desc   = '2D Elliptical Domain - Interstitial';

% --- SCENARIO 2: 3D Compact Box ---
scn2.id = 4;
scn2.type = 'compact';
scn2.dim_mode = '3D';
scn2.shape = 'box';
scn2.angle = [0, 0];
scn2.density = 0.35;
scn2.domain = [0.01, 2.0, 1.0, 0.5]; % [dx, Lx, Ly, Lz]
scn2.core   = [1.5, 0.8, 0.3];       % [width, height, depth]
scn2.desc   = '3D Structured Box (Compact)';

test_scenarios = {scn1, scn2};
seed = 2026;
save_mesh = true;
save_figure = true;

for i = 1:length(test_scenarios)
    s = test_scenarios{i};
    
    fname = sprintf('id%d_%s_%s_%d', s.id, s.dim_mode, s.type, seed);
    
    fprintf('>>> Running Scenario %d: %s\n', i, s.desc);
    cd('..');
    run_fibrosis_generator(s.type, s.density, seed, s.angle, s.dim_mode, ...
                           s.domain, s.shape, s.core, fullfile('examples/outputs', fname), save_mesh, save_figure);
    cd('examples');
end