function [permute_table, offset_table1, offset_table2] = generateTables(seed, is3D)
% GENERATETABLES Pre-computes random permutation arrays and spatial offsets for Perlin Noise.
%
% Generates deterministic mathematical tables based on a user-provided seed. 
% This guarantees absolute reproducibility of the generated histological patterns.
%
% INPUTS:
%   seed - Integer RNG seed.
%   is3D - (Optional) Boolean flag. True generates 3-dimensional offsets.
%
% OUTPUTS:
%   permute_table - [N_freqs x 256] Integer permutation array.
%   offset_table1 - [N_freqs x Dims] Spatial offsets to break grid symmetry.
%   offset_table2 - [N_freqs x Dims] Secondary offsets for advanced flow noise routines.

    if nargin < 2
        is3D = false;
    end

    rng(seed);
    N_freqs = 8;  % Number of octaves/frequencies

    permute_table = zeros(N_freqs, 256, 'int32');
    for j = 1:N_freqs
        permute_table(j,:) = int32(randperm(256) - 1);
    end

    num_cols = 2;
    if is3D
        num_cols = 3;
    end

    offset_table1 = rand(N_freqs, num_cols) - 0.5;
    
    if nargout > 2
        offset_table2 = rand(N_freqs, num_cols) - 0.5;
    else
        offset_table2 = [];
    end
end