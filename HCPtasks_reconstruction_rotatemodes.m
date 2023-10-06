function [cr_original, cr_rotated] = HCPtasks_reconstruction_rotatemodes(X, y, sph, cortex, parcel, n_rep)
% reconstruct after rotating the modes

cr_original = nan(n_rep, 1);
cr_rotated = nan(n_rep, 1);

for i=1:n_rep
    R = rotate_vertices(X, sph.coord, []);
    % set same locations to Zero    
    z = X(:,1)==0 | R(:,1)==0;
    R(z,:) = NaN;
    Xz = X;
    Xz(z,:) = NaN;
    
    [lm, fitted, cr] =  fitlm_cortex(Xz, y, cortex, parcel);
    cr_original(i) = cr.rho;
    [lm, fitted, cr] =  fitlm_cortex(R, y, cortex, parcel);
    cr_rotated(i) = cr.rho;
end
