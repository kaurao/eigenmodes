%% set up
addpath(genpath(fullfile('BrainSpace-0.1.10','matlab')))
addpath(genpath('matlab_GIfTI-master'))

seed = 1;
num_modes = 200;
parcellation_eval = 'Glasser360';

% surface data fron BrainSpace
[surf_lh, ~] = load_conte69();

% get Pang et al.'s data
pang_dir = '../BrainEigenmodes';
parcellations = get_pang_parcellations(pang_dir);
eigenmodes = get_pang_eigenmodes(pang_dir);
eigenmodes_names = fields(eigenmodes);
eigenmodes_l2norm = [];
% filter out modes if needed
for i_emode=1:length(eigenmodes_names)
    emode = eigenmodes_names{i_emode};
    tmp = eigenmodes.(emode);
    eigenmodes.(emode) = tmp(:,1:num_modes);
    eigenmodes_l2norm.(emode) = vecnorm(eigenmodes.(emode),2,2);
    h = plot_hemispheres(eigenmodes_l2norm.(emode), surf_lh);
end

hcptask_file = fullfile(pang_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat');
load(hcptask_file) % this loads `zstat`
all_tasks = fields(zstat);
key_tasks = {'emotion_faces_shapes', 'gambling_punish_reward', ...
            'language_math_story', 'motor_cue_avg', 'relational_match_rel', ...
            'social_tom_random', 'wm_2bk_0bk'};

use_tasks = all_tasks;
zstat_avg  = [];
for i_task=1:length(use_tasks)
    task = use_tasks{i_task};    
    zstat_avg.(task) = mean(zstat.(task), 2, 'omitnan');
end
zstat = []; % save memory
        
eigenmodes_names = fields(eigenmodes);
% original task maps
recon_err  = [];
for i_task=1:length(use_tasks)
    task = use_tasks{i_task};
    fprintf(1, 'task %d/%d: %s\n', i_task, length(use_tasks), task)
    task_avg = zstat_avg.(task);
    recon_err.(task) = zeros(size(eigenmodes.(eigenmodes_names{1}),1), length(eigenmodes_names));
    for i_emode=1:length(eigenmodes_names)
        emode = eigenmodes_names{i_emode};
        [lm, fitted, cr] = fitlm_cortex(eigenmodes.(emode), task_avg, ...
            parcellations.cortex, parcellations.(parcellation_eval));
        recon_err.(task)(:,i_emode) = task_avg - fitted;
    end
end

% concatename across task maps
l2_err = [];
corrs = [];
parcellation = parcellations.(parcellation_eval);
cortex = parcellations.cortex;
for i_emode=1:length(eigenmodes_names)
    emode = eigenmodes_names{i_emode};
    corrs.(emode) = [];
    l2_err.(emode) = [];
    for i_task=1:length(use_tasks)
        task = use_tasks{i_task};
        l2_err_i = [];
        l2_err_i(:,1) = mean_parcel_cortex(eigenmodes_l2norm.(emode), parcellation, cortex);             
        l2_err_i(:,2) = mean_parcel_cortex(recon_err.(task)(:,i_emode), parcellation, cortex);
        l2_err.(emode) = [l2_err.(emode); l2_err_i];
        corrs.(emode)(i_task) = corr(l2_err_i(:,1), abs(recon_err.(task)(:,i_emode)), 'Rows', 'complete');
    end
    cr = corr(l2_err.(emode)(:,1), abs(l2_err.(emode)(:,2)), 'Rows', 'complete');
    fprintf(1, 'Mode: %s, Corr: %.4f\n', emode, cr)
    figure;
    binscatter(l2_err.(emode)(:,1), abs(l2_err.(emode)(:,2)));
    title(emode);
    xlabel('L2-norm'); ylabel('abs(Reconstruction error)');    
end
