% EXAMPLE: PATIENT-SPECIFIC MESH GENERATION
% Demonstrates how to inject procedural fibrosis into custom anatomical models.

clc;
fprintf("=== Cardiac Fibrosis Framework Example: Custom Patient Meshes ===\n\n");

% Configure these paths to point to your local patient datasets
filepath_2D = '../patient_meshes/200um/outputs_2D_dx0.2/Patient_7/Patient_7_slice_6.vtu';
filepath_3D = '../patient_meshes/200um/outputs_3D_dx0.2/Patient_7/Patient_7.alg';

% --- SCENARIO 1: 2D Custom Mesh (VTU format) ---
scn1.type = 'diffuse'; 
scn1.dim_mode = 'CUSTOM'; 
scn1.shape = 'custom'; 
scn1.angle = 45;       
scn1.density = 0.40;   
scn1.domain = filepath_2D; 
scn1.core   = [];      
scn1.desc   = 'Patient (2D Slice) - Diffuse';

% --- SCENARIO 2: 3D Custom Mesh (ALG format) ---
scn2.type = 'patchy'; 
scn2.dim_mode = 'CUSTOM'; 
scn2.shape = 'custom';   
scn2.angle = [0, 0];   
scn2.density = 0.30;   
scn2.domain = filepath_3D; 
scn2.core   = [];      
scn2.desc   = 'Patient (3D Heart) - Patchy';

test_scenarios = {scn1, scn2};
seed = 2026;
save_mesh = true;   
save_figure = false; 

% Execution Loop
for i = 1:length(test_scenarios)
    s = test_scenarios{i};
    fname = sprintf('output_patient_scenario_%d', i);
    
    fprintf('>>> Running Scenario %d: %s\n', i, s.desc);
    % Navigate back to root to call the generator
    cd('..');
    run_fibrosis_generator(s.type, s.density, seed, s.angle, s.dim_mode, ...
                           s.domain, s.shape, s.core, fullfile('examples/outputs', fname), save_mesh, save_figure);
    cd('examples');
end