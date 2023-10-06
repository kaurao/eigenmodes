function [zstat_rotates, zstat_avg] = generate_rotates(randomization_with, randomization_on, num_rotates, seed, ...
    pang_dir, key_tasks, result_dir) 
% randomization_on: 'spheres' or 'surfaces', used when randomization_with is 'spin'
% randomization_with: 'spin' or 'moran'
% needs BrainSpace and matlab_GIfTI in path

%% set up
%addpath(genpath(fullfile('..', 'BrainSpace','matlab')))
%addpath(genpath(fullfile('..', 'matlab_GIfTI')))
if ~exist('spin_permutations', 'file')
    error('spin_permutations not found, please add BrainSpace to path')
end

if ~exist('pang_dir', 'var') || isempty(pang_dir)
    pang_dir = fullfile('..','BrainEigenmodes');
end
parcellations = get_pang_parcellations(pang_dir);

if ~exist('result_dir', 'var')
    result_dir = '';
end
if ~isempty(result_dir) && ~exist(result_dir, 'dir')
    mkdir(result_dir)
elseif isempty(result_dir)
    fprintf(1, 'result_dir not specified, will not save results\n')
end

if ~exist('seed', 'var') || isempty(seed)
    seed = 1;
end
if ~exist('num_rotates', 'var') || isempty(num_rotates)
    num_rotates = 5000;
end

if ~exist('key_tasks', 'var') || isempty(key_tasks)
    key_tasks = {'emotion_faces_shapes', 'gambling_punish_reward', ...
            'language_math_story', 'motor_cue_avg', 'relational_match_rel', ...
            'social_tom_random', 'wm_2bk_0bk'};
end

if strcmp(lower(randomization_with), 'spin')
    assert(exist('randomization_on', 'var') && ~isempty(randomization_on), ...
        'randomization_on must be specified when randomization_with is spin')
end

fprintf(1, 'generating %d rotations: %s %s\n', ...
    num_rotates, randomization_on, randomization_with)

%% load data and compute rotations
hcptask_file = fullfile(pang_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat');
load(hcptask_file) % this loads `zstat`

zstat_avg  = [];
for i_task=1:length(key_tasks)
    task = key_tasks{i_task};    
    zstat_avg.(task) = mean(zstat.(task), 2, 'omitnan');
end
zstat = []; % save memory

zstat_rotates = [];
[sphere_lh, ~] = load_conte69(randomization_on);
if lower(randomization_with) == 'moran'
    filename = sprintf('MEM_%s.mat', randomization_on);
    filename = fullfile(result_dir, filename);
    if exist(filename)
        fprintf(1, 'loading %s\n', filename)
        load(filename)
    else
        fprintf(1, 'running compute_mem for moran... ')
        % we are using parameters as used by Faskowitz et al
        MEM = compute_mem(sphere_lh, 'mask', ~parcellations.cortex, ...
        'n_ring', 1e4, 'distance_metric','geodesic');
        save(filename, 'MEM', '-v7.3')
        fprintf('done\n')
    end
end

for i_task=1:length(key_tasks)
    tic
    task = key_tasks{i_task};
    fprintf(1, 'task %d/%d: %s\n', i_task, length(key_tasks), task)
    task_avg = zstat_avg.(task);
    rng(seed) % keep it reproducible
    if lower(randomization_with) == 'moran'
        task_rotates = moran_randomization(task_avg(parcellations.cortex),...
                                MEM,num_rotates,'procedure','singleton');
    elseif lower(randomization_with) == 'spin'
        task_rotates = spin_permutations(task_avg,sphere_lh,num_rotates);
    else
        error(['unknown randomization: ' randomization_with])
    end
    zstat_rotates.(task) = task_rotates;
    toc
end


if ~isempty(result_dir)
    filename = sprintf('HCPtasks_%drotates_rotation_%s-%s_seed%d.mat', ...
        num_rotates, randomization_on, randomization_with, seed);
    result_file = fullfile(result_dir, filename);
    save(result_file, 'seed', 'key_tasks', ...
        'zstat_avg', 'zstat_rotates', '-v7.3')
end
