function [lm, fitted, cr] =  fitlm_cortex(X, y, cortex, parcellation)
% fit linear model to cortex vertices
% the function takes care of nan entries in y
% X: feature matrix, e.g. eigenmodes
% y: target, e.g. task data
% cortex: logical vector with 1 for cortex vertices
% parcellation: parcellation scheme for correlation
% written by Kaustubh R. Patil, 2023

cr = [];
cr.rho = nan;
cr.pval = nan;

if ~exist('parcellation', 'var')
    parcellation = [];
end

% take care of possible nan entries in y
% for instance they can be spun medial wall
cortex = cortex & ~isnan(y);

lm = fitlm(X(cortex,:), y(cortex));
fitted = nan(length(y), 1);
fitted(cortex) = lm.Fitted;

% get correlation based on provided parcellation scheme
if ~isempty(parcellation)
    yparc = mean_parcel_cortex(y, parcellation, cortex);
    fparc = mean_parcel_cortex(fitted, parcellation, cortex);
    [cr.rho, cr.pval] = corr(yparc, fparc, 'Rows','complete');
end

