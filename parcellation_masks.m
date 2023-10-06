function [parc_masks, u_parcels] = parcellation_masks(parcellation)
% given a parcellation create basis vectors from it
% parcellation: a vector with each element corresponding to a vertex
% parc_masks: a matrix with each column corresponding to a parcel
% u_parcels: unique parcels
% written by Kaustubh R. Patil, 2023

u_parcels = unique(parcellation);
u_parcels(u_parcels<=0) = [];
u_parcels(isnan(u_parcels)) = [];
parc_masks = zeros(length(parcellation),length(u_parcels));
for i=1:length(u_parcels)
    idx = parcellation==u_parcels(i);
    parc_masks(idx, i) = 1;
end
