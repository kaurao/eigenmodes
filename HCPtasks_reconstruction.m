addpath(genpath(fullfile('BrainSpace-0.1.10','matlab')))
addpath(genpath('matlab_GIfTI-master'))

if ~exist('num_modes','var')
    num_modes = 200;
end

parcellation_mode = 'Schaefer';  % 'Schaefer' or 'kmedoids'
parcellation_eval = 'Glasser360';

% get Pang et al.'s data
pang_dir = '../BrainEigenmodes';
parcellations = get_pang_parcellations(pang_dir);
eigenmodes = get_pang_eigenmodes(pang_dir);
eigenmodes_names = fields(eigenmodes);
for i_emode=1:length(eigenmodes_names)
    emode = eigenmodes_names{i_emode};
    tmp = eigenmodes.(emode);
    eigenmodes.(emode) = tmp(:,1:num_modes);
end

%load `avgSC_L`
load(fullfile(pang_dir, 'data', 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'));
avgSC_L = log(avgSC_L+1);

hcptask_file = fullfile(pang_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat');
load(hcptask_file) % this loads `zstat`
key_tasks = {'emotion_faces_shapes', 'gambling_punish_reward', ...
            'language_math_story', 'motor_cue_avg', 'relational_match_rel', ...
            'social_tom_random', 'wm_2bk_0bk'};
all_tasks = fields(zstat);

zstat_avg  = [];
for i_task=1:length(all_tasks)
    task = all_tasks{i_task};    
    zstat_avg.(task) = mean(zstat.(task), 2, 'omitnan');
end
zstat = []; % save memory
        
% sufcae data fron BrainSpace
[surf_lh, ~] = load_conte69();

% generate parcel-informed basis vectors
if strcmp(parcellation_mode, 'Schaefer')
    parcellation_mode = sprintf('Schaefer%d', num_modes*2);
elseif strcmp(parcellation_mode, 'kmedoids')
    parcellation_mode = sprintf('kmedoids%d', num_modes);
    parc = readmatrix(sprintf('../kmedoids%d_geodesic_lh.txt', num_modes));
    parc(isnan(parc)) = 0;
    parcellations.(parcellation_mode) = parc;
end
[eigenmodes.Schaefer_mask, u_parcels] = parcellation_masks(parcellations.(parcellation_mode));
eigenmodes.Schaefer_Connectome = eigenmodes.Schaefer_mask;
eigenmodes.Schaefer_Connectome(:) = 0;
parc = parcellations.(parcellation_mode);
parc_ctx = parc(parcellations.cortex);
for i=1:size(eigenmodes.Schaefer_mask,2)
   eigenmodes.Schaefer_Connectome(parcellations.cortex,i) = mean(avgSC_L(parc_ctx==i,:), 'omitnan')';
   eigenmodes.Schaefer_Connectome(parc==i,i) = max(eigenmodes.Schaefer_Connectome(parc==i,i));
end

eigenmodes_names = fields(eigenmodes);
recon_acc = nan(length(all_tasks), length(eigenmodes_names));
for i_task=1:length(all_tasks)
    task = all_tasks{i_task};
    fprintf(1, 'task %d/%d: %s\n', i_task, length(all_tasks), task)    
    task_avg = zstat_avg.(task);
    lm = [];
    fitted = [];
    cr  = [];
    for i_emode=1:length(eigenmodes_names)
        emode = eigenmodes_names{i_emode};
        [lm.(emode), fitted.(emode), cr.(emode)] = fitlm_cortex(eigenmodes.(emode), task_avg, ...
                        parcellations.cortex, parcellations.(parcellation_eval));
        recon_acc(i_task, i_emode) = cr.(emode).rho;
    end
end

% boxplot(recon_acc, 'Labels', eigenmodes_names)
