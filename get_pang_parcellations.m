function parcellations_lh = get_pang_parcellations(pang_dir)
% get parcellations as provided by Pang et al. 2023
% returns a structure with Schaefer parcellation and cortex indices
% works only for lh
% written by Kaustubh R. Patil, 2023

parcellation_dir = fullfile(pang_dir, 'data', 'parcellations');
parcellation_files = dir(fullfile(parcellation_dir, '*-lh.txt'));

parcellations_lh = [];
for i=1:length(parcellation_files)
    name = parcellation_files(i).name;
    [~, name_short, ~] = fileparts(name);
    name_short = strsplit(name_short, '_');
    name_short = name_short{end};
    name_short = strsplit(name_short, '-');
    name_short = name_short{1};
    parcellations_lh.(name_short) = int16(dlmread(fullfile(parcellation_dir, name)));
end

% get cortex indices
parcellations_lh.cortex = dlmread(fullfile(pang_dir, 'data', 'template_surfaces_volumes', 'fsLR_32k_cortex-lh_mask.txt'));
parcellations_lh.cortex = logical(parcellations_lh.cortex);
