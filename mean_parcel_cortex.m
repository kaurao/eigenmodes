function xmean = mean_parcel_cortex(x, parcellation, cortex)
% get mean of each parcel
% x: data vector
% parcellation: parcellation vector
% cortex: logical vector with 1 for cortex vertices
	assert(length(parcellation)==length(cortex))
	u_parcels = unique(parcellation(cortex));
	u_parcels(u_parcels<=0) = [];
    xmean = nan(length(u_parcels),1);
	for i=1:length(u_parcels)
		idx = parcellation==u_parcels(i);
		xmean(i) = mean(x(idx), 'omitnan');
	end

