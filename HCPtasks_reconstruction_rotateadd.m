% this script reconstructs the seven HCP key tasks using specified number of modes
% and compares the results with the original task maps
% expects the data from Pang et al. (2023) to be in `pang_dir`
% also BrainSpace and matlab_GIfTI should be in the path
% parallel computing toolbox is needed for the parfor loop
% written by Kaustubh R. Patil, 2023

%% set up
addpath(genpath(fullfile('..', 'BrainSpace','matlab')))
addpath(genpath(fullfile('..', 'matlab_GIfTI')))

pang_dir = '../BrainEigenmodes';

seed = 1; % for reproducibility
num_modes = 200; % number of modes to use
num_rotates = 10; % number of rotations to use, 10 is for demonstration use 5000 for the paper
parcellation_mode = 'Schaefer'; % 'Schaefer' or 'kmedoids'
parcellation_eval = 'Glasser360'; % for evaluation
randomization_on = 'spheres'; % used for 'spin': 'spheres' or 'surfaces' ('spheres' is standard)

% Uncomment these two lines and set values as required to use this as a standalone script
% ---------------------------------------------------------------------------------------
% randomization_with = 'moran'; % 'moran' or 'spin' or 'file'
% randomization_file = '../results/HCPtasks_5000rotates_rotation_surfaces-moran_seed1.mat';
if strcmp(lower(randomization_with), 'file')
    fprintf(1, 'loading rotated data: %s\n', randomization_file)
    assert(exist(randomization_file, 'file'))
    load(randomization_file);
end

% get Pang et al. data
parcellations = get_pang_parcellations(pang_dir);
eigenmodes = get_pang_eigenmodes(pang_dir);
eigenmodes_names = fields(eigenmodes);
% filter out modes if needed
for i_mode=1:length(eigenmodes_names)
    emode = eigenmodes_names{i_mode};
    tmp = eigenmodes.(emode);
    eigenmodes.(emode) = tmp(:,1:num_modes);
end

hcptask_file = fullfile(pang_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat');
load(hcptask_file) % this loads `zstat`
key_tasks = {'emotion_faces_shapes', 'gambling_punish_reward', ...
            'language_math_story', 'motor_cue_avg', 'relational_match_rel', ...
            'social_tom_random', 'wm_2bk_0bk'};

% get average task maps
zstat_avg  = [];
for i_task=1:length(key_tasks)
    task = key_tasks{i_task};    
    zstat_avg.(task) = mean(zstat.(task), 2, 'omitnan');
end
zstat = []; % save memory
        
%% parcel-informed basis vectors
% corresponding Schaefer parcels 
% generate parcel-informed basis vectors
if strcmp(parcellation_mode, 'Schaefer')
    parcellation_mode = sprintf('Schaefer%d', num_modes*2);
elseif strcmp(parcellation_mode, 'kmedoids')
    parcellation_mode = sprintf('kmedoids%d', num_modes);
    parc = readmatrix(sprintf('../kmedoids%d_geodesic_lh.txt', num_modes));
    parc(isnan(parc)) = 0;
    parcellations.(parcellation_mode) = parc;
else
    error(['unknown parcellation_mode: ' parcellation_mode])
end

% get Parcel basis vectors
[eigenmodes.Schaefer_mask, u_parcels] = parcellation_masks(parcellations.(parcellation_mode));
% generate Parcel+Connectome basis vectors by adding SC to parcels
load(fullfile(pang_dir, 'data', 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'));
avgSC_L = log(avgSC_L+1);
eigenmodes.Schaefer_Connectome = eigenmodes.Schaefer_mask;
eigenmodes.Schaefer_Connectome(:) = 0;
parc = parcellations.(parcellation_mode);
parc_ctx = parc(parcellations.cortex);
for i=1:size(eigenmodes.Schaefer_mask,2)
   eigenmodes.Schaefer_Connectome(parcellations.cortex,i) = mean(avgSC_L(parc_ctx==i,:), 'omitnan')';
   eigenmodes.Schaefer_Connectome(parc==i,i) = max(eigenmodes.Schaefer_Connectome(parc==i,i));
end
avgSC_L = []; % save memory
eigenmodes_names = fields(eigenmodes);

%% reconstruct original task maps
recon_orig = zeros(length(key_tasks), length(eigenmodes_names));
recon_err  = zeros(length(key_tasks), size(eigenmodes.(eigenmodes_names{1}),1));
for i_task=1:length(key_tasks)
    task = key_tasks{i_task};
    fprintf(1, 'task %d/%d: %s\n', i_task, length(key_tasks), task)
    task_avg = zstat_avg.(task);
    for i_mode=1:length(eigenmodes_names)
        emode = eigenmodes_names{i_mode};
        [lm, fitted, cr] = fitlm_cortex(eigenmodes.(emode), task_avg, ...
            parcellations.cortex, parcellations.(parcellation_eval));
        recon_orig(i_task, i_mode) = cr.rho;
    end
end

%% reconstruct randomzied task maps
recon_rotateadd = [];
recon_rotate = [];
corr_geom_beta = [];
corr_original_rotate = [];
[sphere_lh, ~] = load_conte69(randomization_on);

if strcmp(lower(randomization_with), 'moran')
    fprintf(1, 'running compute_mem for moran... ')
    MEM = compute_mem(sphere_lh, 'mask', ~parcellations.cortex);
        %'n_ring', 1e4, 'distance_metric','mesh');
    fprintf('done\n')
end
num_modes = length(eigenmodes_names);
for i_task=1:length(key_tasks)
    recon_r = nan(num_rotates, num_modes);
    recon_ra = nan(num_rotates, num_modes);
    task = key_tasks{i_task};
    fprintf(1, 'task %d/%d: %s\n', i_task, length(key_tasks), task)
    task_avg = zstat_avg.(task);
    lm = fitlm_cortex(eigenmodes.Geometric, task_avg, ...
                        parcellations.cortex, parcellations.(parcellation_eval));
    coeff = lm.Coefficients.Estimate(3:end); % 1 is Intercept and also 2
    corr_gb = zeros(1, num_rotates); % beta
    corr_or = zeros(1, num_rotates); % map
    rng(seed) % keep it reproducible
    if strcmp(lower(randomization_with), 'moran')
        task_rotates = moran_randomization(task_avg(parcellations.cortex),...
                                MEM,num_rotates,'procedure','singleton');
    elseif strcmp(lower(randomization_with), 'spin')
        task_rotates = spin_permutations(task_avg,sphere_lh,num_rotates);
    elseif strcmp(lower(randomization_with), 'file')
        task_rotates = zstat_rotates.(task);
    else
        error(['unknown randomization: ' randomization_with])
    end

    % take care of different data structures generated by rotations
    if iscell(task_rotates)
        assert(numel(task_rotates)==1)
        task_rotates = task_rotates{1};
    end

    parfor i_rot=1:num_rotates        
        task_rand = task_rotates(:,:,i_rot);
        % make sure that the rotated tasks have all vertices
        if size(task_rand,1)==29696
            tmp_task_rand = nan(size(task_avg));
            tmp_task_rand(parcellations.cortex) = task_rand;
            task_rand = tmp_task_rand;
        end
        corr_or(i_rot) = corr(task_avg, task_rand, 'Rows', 'complete');
        
        % to follow Faskowitz et al. we will remove nan values in
        % fitlm_cortex and correlate with 'Rows','complete'
        %task_rand(isnan(task_rand)) = 0;
        task_avg_rand = task_avg + task_rand;        
        for i_mode=1:num_modes
            emode = eigenmodes_names{i_mode};
            [lm_ra, fitted_ra, cr_ra] = fitlm_cortex(eigenmodes.(emode), task_avg_rand, ...
                        parcellations.cortex, parcellations.(parcellation_eval));
            recon_ra(i_rot, i_mode) = cr_ra.rho;
            
            [lm_r, fitted_r, cr_r] = fitlm_cortex(eigenmodes.(emode), task_rand, ...
                        parcellations.cortex, parcellations.(parcellation_eval));
            recon_r(i_rot, i_mode) = cr_r.rho;
            
            % get the correlation between Geometric beta values
            if strmp(lower(emode), 'geometric')
                coeff_rand = lm_r.Coefficients.Estimate(3:end); % remove 1 and 2
                corr_gb(i_rot) = corr(coeff, coeff_rand, 'Rows', 'Complete');
            end
        end
    end
    recon_rotateadd.(task) = recon_ra;
    recon_rotate.(task) = recon_r;
    corr_geom_beta.(task) = corr_gb;
    corr_original_rotate.(task) = corr_or;
end
