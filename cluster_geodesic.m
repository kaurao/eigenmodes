function clusters = cluster_geodesic(geodesic, cortex, K)
% this function clusters the vertices of the fsLR 32k surface based on geodesic distance
% the nan entries depicting the medial wall is excluded before clustering
% geodesic: geodesic distance matrix
%           e.g. load fsLR_32k_high-resolution_geodesic_dist_midthickness-lh.mat
% cortex: logical vector with 1 for cortex vertices
%         e.g. cortex = readmatrix('fsLR_32k_cortex-lh_mask.txt');
% K: number of clusters
% written by Kaustubh R. Patil, 2023

% remove medial wall
geodesic = geodesic(cortex>0, cortex>0);
% make it symmetric
geodesic = geodesic + geodesic';

N = size(geodesic,1);
% as kmedoids cannot handle precomputed distance matrix
% we use a trick to provide it with a distance function
[idx, C, sumd] = kmedoids((1:N)', K, 'Distance', @(ZI, ZJ) geodesic(ZJ, ZI));

clusters = nan(length(cortex),1);
clusters(cortex>0) = idx;