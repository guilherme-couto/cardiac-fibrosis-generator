% EXAMPLE: PATIENT-SPECIFIC MESH GENERATION
% Demonstrates how to inject procedural fibrosis into custom anatomical models.

clc;
fprintf("Cardiac Fibrosis Framework Example: Custom Meshes\n\n");

% Configure these paths to point to your local files
filepath_2D = '/path/to/your/patient_mesh_2D_slice.alg';  % Replace with your actual 2D slice mesh file path
filepath_3D = '/path/to/your/patient_mesh_3D.alg';  % Replace with your actual 3D mesh file path

% --- SCENARIO 1: 2D Custom Mesh ---
scn1.id = 1;
scn1.type = 'diffuse'; 
scn1.dim_mode = 'CUSTOM'; 
scn1.shape = []; 
scn1.angle = [];       
scn1.density = 0.40;   
scn1.domain = filepath_2D; 
scn1.core   = [];      
scn1.desc   = 'Patient (2D Slice) - Diffuse';

% --- SCENARIO 2: 3D Custom Mesh ---
scn2.id = 2;
scn2.type = 'patchy'; 
scn2.dim_mode = 'CUSTOM'; 
scn2.shape = [];   
scn2.angle = []; % Will be ignored  
scn2.density = 0.30;   
scn2.domain = filepath_3D; 
scn2.core   = [];      
scn2.desc   = 'Patient (3D Heart) - Patchy';

test_scenarios = {scn1, scn2};
seed = 2026;
save_mesh = true;   
save_figure = true; 

% Execution Loop
for i = 1:length(test_scenarios)
    s = test_scenarios{i};

    fname = sprintf('id%d_%s_%s_%d', s.id, s.dim_mode, s.type, seed);
    
    fprintf('>>> Running Scenario %d: %s\n', i, s.desc);
    % Navigate back to root to call the generator
    cd('..');
    run_fibrosis_generator(s.type, s.density, seed, s.angle, s.dim_mode, ...
                           s.domain, s.shape, s.core, fullfile('examples/outputs', fname), save_mesh, save_figure);
    cd('examples');
end