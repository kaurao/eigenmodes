function eigenmodes_lh = get_pang_eigenmodes(pang_dir)
% get eigenmodes as provided by Pang et al. 2023
% returns a structure with fields:
%   Geometric: geometric eigenmodes
%   Connectome: connectome eigenmodes
%   Connectome_DM: connectome eigenmodes with density matching
%   EDR: eigenmodes from EDR
% written by Kaustubh R. Patil, 2023

emode_dir = fullfile(pang_dir, 'data', 'results');

eigenmodes_lh = [];

eigenmodes_lh.Geometric = dlmread(fullfile(emode_dir, 'basis_geometric_midthickness-lh_evec_200.txt'));
temp = load(fullfile(emode_dir, 'basis_connectome_midthickness-lh_evec_200.mat'));
eigenmodes_lh.Connectome = temp.eig_vec;
temp = load(fullfile(emode_dir, 'basis_connectome_density_matched_midthickness-lh_evec_200.mat'));
eigenmodes_lh.Connectome_DM = temp.eig_vec;
temp = load(fullfile(emode_dir, 'basis_connectome_EDR_midthickness-lh_evec_200.mat'));
eigenmodes_lh.EDR = temp.eig_vec;
