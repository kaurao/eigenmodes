% Geodesic distance
% -----------------
if ~isfile('./results/fsLR_32k_high-resolution_geodesic_dist_midthickness-lh.mat')
    hemisphere = 'lh';
    mesh_interest = 'midthickness';
    [vertices, faces] = read_vtk(sprintf('../BrainEigenmodes/data/template_surfaces_volumes/fsLR_32k_%s-%s.vtk', mesh_interest, hemisphere));
    surface_midthickness.vertices = vertices';
    surface_midthickness.faces = faces';
    gSurface = gifti(surface_midthickness);
    save(gSurface,'./results/surface_midthickness_lh.gii');
    
    % Workbench 
    % ------------------------------
    % High resolution 32492 vertices
    % ------------------------------
    unix('wb_command -surface-geodesic-distance-all-to-all ./results/surface_midthickness_lh.gii ./results/geodesic_dist_midthickness_lh.dconn.nii');
    unix('wb_command -cifti-convert -to-text ./results/geodesic_dist_midthickness_lh.dconn.nii ./results/geodesic_dist_midthickness_lh.dconn.csv -col-delim ,');
    geodesic_lh = dlmread('./results/geodesic_dist_midthickness_lh.dconn.csv');
    save('./results/fsLR_32k_high-resolution_geodesic_dist_midthickness-lh.mat','geodesic_lh');
end

% kmedoids
% --------
if ~isfile('./results/kmedoids50_geodesic_lh.txt') || ~isfile('./results/kmedoids100_geodesic_lh.txt') || ~isfile('./results/kmedoids180_geodesic_lh.txt') || ~isfile('./results/kmedoids200_geodesic_lh.txt')
    load('./results/fsLR_32k_high-resolution_geodesic_dist_midthickness-lh.mat','geodesic_lh')
    rng('default');
    geodesic_lh = geodesic_lh(cortex>0, cortex>0);
    geodesic_lh = geodesic_lh + geodesic_lh';
    N = size(geodesic_lh,1);
end
if isfile('./results/kmedoids50_geodesic_lh.txt')
    atlasfree50  = dlmread('./results/kmedoids50_geodesic_lh.txt');
else
    K = 50;% number of clusters
    [idx, C, sumd] = kmedoids((1:N)', K, 'Distance', @(ZI, ZJ) geodesic_lh(ZJ, ZI));
    clusters = nan(length(cortex),1);
    clusters(cortex>0) = idx;
    dlmwrite('./results/kmedoids50_geodesic_lh.txt',clusters);
    atlasfree50 = clusters;
end
if isfile('./results/kmedoids100_geodesic_lh.txt')
    atlasfree100  = dlmread('./results/kmedoids100_geodesic_lh.txt');
else
    K = 100;
    [idx, C, sumd] = kmedoids((1:N)', K, 'Distance', @(ZI, ZJ) geodesic_lh(ZJ, ZI));
    clusters = nan(length(cortex),1);
    clusters(cortex>0) = idx;
    dlmwrite('./results/kmedoids100_geodesic_lh.txt',clusters);
    atlasfree100 = clusters;
end
if isfile('./results/kmedoids180_geodesic_lh.txt')
    atlasfree180  = dlmread('./results/kmedoids180_geodesic_lh.txt');
else
    K = 180;
    [idx, C, sumd] = kmedoids((1:N)', K, 'Distance', @(ZI, ZJ) geodesic_lh(ZJ, ZI));
    clusters = nan(length(cortex),1);
    clusters(cortex>0) = idx;
    dlmwrite('./results/kmedoids180_geodesic_lh.txt',clusters);
    atlasfree180 = clusters;
end
if isfile('./results/kmedoids200_geodesic_lh.txt')
    atlasfree200  = dlmread('./results/kmedoids200_geodesic_lh.txt');
else
    K = 200;
    [idx, C, sumd] = kmedoids((1:N)', K, 'Distance', @(ZI, ZJ) geodesic_lh(ZJ, ZI));
    clusters = nan(length(cortex),1);
    clusters(cortex>0) = idx;
    dlmwrite('./results/kmedoids200_geodesic_lh.txt',clusters);
    atlasfree200 = clusters;
end

% Reconstruction accuracies (parcellations)
% -----------------------------------------
if isfile('./results/reconstruction_accuracy.mat')
    load('./results/reconstruction_accuracy.mat','eigenmodes_n50','eigenmodes_n100','eigenmodes_n180','eigenmodes_n200','recon_acc_n50','recon_acc_n100','recon_acc_n180','recon_acc_n200','eigenmodes_names');
else
    num_modes = 200;
    parcellation_mode = 'Schaefer';
    HCPtasks_reconstruction;
    eigenmodes_n200  = eigenmodes;
    recon_acc_n200   = recon_acc;
    eigenmodes_names = eigenmodes_names;
    
    num_modes = 180;
    parcellation_mode = 'Glasser';
    HCPtasks_reconstruction;
    eigenmodes_n180  = eigenmodes;
    recon_acc_n180   = recon_acc;
    eigenmodes_names = eigenmodes_names;
    
    num_modes = 100;
    parcellation_mode = 'Schaefer';
    HCPtasks_reconstruction;
    eigenmodes_n100  = eigenmodes;
    recon_acc_n100   = recon_acc;
    eigenmodes_names = eigenmodes_names;
    
    num_modes = 50;
    parcellation_mode = 'Schaefer';
    HCPtasks_reconstruction;
    eigenmodes_n50  = eigenmodes;
    recon_acc_n50   = recon_acc;
    eigenmodes_names = eigenmodes_names;
    save('./results/reconstruction_accuracy.mat','eigenmodes_n200','eigenmodes_n180','eigenmodes_n100','eigenmodes_n50','recon_acc_n200','recon_acc_n180','recon_acc_n100','recon_acc_n50','eigenmodes_names');
end

% Reconstruction accuracies (kmedoids)
% ------------------------------------
if isfile('./results/reconstruction_accuracy_kmedoids.mat')
    load('./results/reconstruction_accuracy_kmedoids.mat','eigenmodes_kmedoids_n200','eigenmodes_kmedoids_n180','eigenmodes_kmedoids_n100','eigenmodes_kmedoids_n50','recon_acc_kmedoids_n200','recon_acc_kmedoids_n180','recon_acc_kmedoids_n100','recon_acc_kmedoids_n50','eigenmodes_names_kmedoids');
else
    num_modes = 200;
    parcellation_mode = 'kmedoids';
    HCPtasks_reconstruction;
    eigenmodes_kmedoids_n200  = eigenmodes;
    recon_acc_kmedoids_n200   = recon_acc;
    eigenmodes_names_kmedoids = eigenmodes_names;
    
    num_modes = 180;
    parcellation_mode = 'kmedoids';
    HCPtasks_reconstruction;
    eigenmodes_kmedoids_n180  = eigenmodes;
    recon_acc_kmedoids_n180   = recon_acc;
    eigenmodes_names_kmedoids = eigenmodes_names;
    
    num_modes = 100;
    parcellation_mode = 'kmedoids';
    HCPtasks_reconstruction;
    eigenmodes_kmedoids_n100  = eigenmodes;
    recon_acc_kmedoids_n100   = recon_acc;
    eigenmodes_names_kmedoids = eigenmodes_names;
    
    num_modes = 50;
    parcellation_mode = 'kmedoids';
    HCPtasks_reconstruction;
    eigenmodes_kmedoids_n50  = eigenmodes;
    recon_acc_kmedoids_n50   = recon_acc;
    eigenmodes_names_kmedoids = eigenmodes_names;
    save('./results/reconstruction_accuracy_kmedoids.mat','eigenmodes_kmedoids_n200','eigenmodes_kmedoids_n180','eigenmodes_kmedoids_n100','eigenmodes_kmedoids_n50','recon_acc_kmedoids_n200','recon_acc_kmedoids_n180','recon_acc_kmedoids_n100','recon_acc_kmedoids_n50','eigenmodes_names_kmedoids');
end

% L2-norm
% -------
if isfile('./results/reconstruction_error.mat')
    load('./results/reconstruction_error.mat','eigenmodes_l2norm','l2_err');
else
    HCPtasks_reconstruction_error;
    save('./results/reconstruction_error.mat','eigenmodes_l2norm','l2_err');
    close all
end

% Spin test
% ---------
if isfile('./results/HCPtasks_5000rotates_seed1.mat')
    load('./results/HCPtasks_5000rotates_seed1.mat','zstat_rotates');
else
    key_tasks = {'emotion_faces_shapes', 'gambling_punish_reward', ...
                 'language_math_story', 'motor_cue_avg', 'relational_match_rel', ...
                 'social_tom_random', 'wm_2bk_0bk'};
    result_dir = fullfile('.', 'results');
    pang_dir = '../BrainEigenmodes';
    [zstat_rotates, zstat_avg] = generate_rotates('spin','spheres',5000,1,pang_dir,key_tasks,result_dir);
    save('./results/HCPtasks_5000rotates_seed1.mat','zstat_rotates','zstat_avg');
end
if isfile('./results/HCPtasks_reconstruction_200modes_5000rotates.mat')
    temp = load('./results/HCPtasks_reconstruction_200modes_5000rotates.mat');
    recon_acc_n200_rot5000 = temp.recon_rotateadd;
    recon_acc_n200_rot5000_no_add = temp.recon_rotate;
    corr_geom_beta_spin = temp.corr_geom_beta;
    corr_orig_rotate_spin = temp.corr_orig_rotate;
    zstat_avg_spin = temp.zstat_avg;
    clear temp
else
    randomization_file = './results/HCPtasks_5000rotates_seed1.mat';
    randomization_with = 'file';
    HCPtasks_reconstruction_rotateadd;
    save('./results/HCPtasks_reconstruction_200modes_5000rotates.mat','recon_rotateadd','recon_rotate','corr_geom_beta','corr_orig_rotate','zstat_avg');
    recon_acc_n200_rot5000 = recon_rotateadd;
    recon_acc_n200_rot5000_no_add = recon_rotate;
    corr_geom_beta_spin = corr_geom_beta;
    corr_orig_rotate_spin = corr_orig_rotate;
    zstat_avg_spin = zstat_avg;
end

% Adjusted Rand index
% -------------------
schaefer400  = dlmread('../BrainEigenmodes/data/parcellations/fsLR_32k_Schaefer400-lh.txt');
schaefer200  = dlmread('../BrainEigenmodes/data/parcellations/fsLR_32k_Schaefer200-lh.txt');
schaefer100  = dlmread('../BrainEigenmodes/data/parcellations/fsLR_32k_Schaefer100-lh.txt');
glasser360   = dlmread('../BrainEigenmodes/data/parcellations/fsLR_32k_Glasser360-lh.txt');
ari200 = rand_index(schaefer400(cortex_ind),atlasfree200(cortex_ind),'adjusted');
ari180 = rand_index(glasser360(cortex_ind)-180,atlasfree180(cortex_ind),'adjusted');
ari100 = rand_index(schaefer200(cortex_ind),atlasfree100(cortex_ind),'adjusted');
ari50  = rand_index(schaefer100(cortex_ind),atlasfree50(cortex_ind),'adjusted');
ari_gla_sch200 = rand_index(glasser360(cortex_ind),schaefer400(cortex_ind),'adjusted');
ari_gla_sch100 = rand_index(glasser360(cortex_ind),schaefer200(cortex_ind),'adjusted');
ari_gla_sch50  = rand_index(glasser360(cortex_ind),schaefer100(cortex_ind),'adjusted');
ari_gla_kmd200 = rand_index(glasser360(cortex_ind),atlasfree200(cortex_ind),'adjusted');
ari_gla_kmd100 = rand_index(glasser360(cortex_ind),atlasfree100(cortex_ind),'adjusted');
ari_gla_kmd50  = rand_index(glasser360(cortex_ind),atlasfree50(cortex_ind),'adjusted');

% Randomized maps (moran)
% -----------------------
if ~isfile('./results/HCPtasks_5000rotates_rotation_surfaces-moran_seed1.mat')
    key_tasks = {'emotion_faces_shapes', 'gambling_punish_reward', ...
                 'language_math_story', 'motor_cue_avg', 'relational_match_rel', ...
                 'social_tom_random', 'wm_2bk_0bk'};
    result_dir = fullfile('.', 'results');
    pang_dir = '../BrainEigenmodes';
    [zstat_rotates, zstat_avg] = generate_rotates('moran',[],5000,1,pang_dir,key_tasks,result_dir);
    save('./results/HCPtasks_5000rotates_rotation_surfaces-moran_seed1.mat','zstat_rotates','zstat_avg');
end
temp = load('./results/HCPtasks_5000rotates_rotation_surfaces-moran_seed1.mat');
l = size(temp.zstat_avg.emotion_faces_shapes,1);
zstat_random.emotion_faces_shapes{1,1}                     = zeros(l,1,5000);
zstat_random.emotion_faces_shapes{1,1}(cortex_ind,1,:)     = temp.zstat_rotates.emotion_faces_shapes;
zstat_random.gambling_punish_reward{1,1}                   = zeros(l,1,5000);
zstat_random.gambling_punish_reward{1,1}(cortex_ind,1,:)   = temp.zstat_rotates.gambling_punish_reward;
zstat_random.language_math_story{1,1}                      = zeros(l,1,5000);
zstat_random.language_math_story{1,1}(cortex_ind,1,:)      = temp.zstat_rotates.language_math_story;
zstat_random.motor_cue_avg{1,1}                            = zeros(l,1,5000);
zstat_random.motor_cue_avg{1,1}(cortex_ind,1,:)            = temp.zstat_rotates.motor_cue_avg;
zstat_random.relational_match_rel{1,1}                     = zeros(l,1,5000);
zstat_random.relational_match_rel{1,1}(cortex_ind,1,:)     = temp.zstat_rotates.relational_match_rel;
zstat_random.social_tom_random{1,1}                        = zeros(l,1,5000);
zstat_random.social_tom_random{1,1}(cortex_ind,1,:)        = temp.zstat_rotates.social_tom_random;
zstat_random.wm_2bk_0bk{1,1}                               = zeros(l,1,5000);
zstat_random.wm_2bk_0bk{1,1}(cortex_ind,1,:)               = temp.zstat_rotates.wm_2bk_0bk;
clear temp
if isfile('./results/HCPtasks_reconstruction_200modes_5000rotates_moran.mat')
    temp = load('./results/HCPtasks_reconstruction_200modes_5000rotates_moran.mat');
    recon_acc_n200_rnd5000 = temp.recon_rotateadd;
    recon_acc_n200_rnd5000_no_add = temp.recon_rotate;
    corr_geom_beta_moran = temp.corr_geom_beta;
    corr_orig_rotate_moran = temp.corr_orig_rotate;
    zstat_avg_moran = temp.zstat_avg;
    clear temp
else
    randomization_file = './results/HCPtasks_5000rotates_rotation_surfaces-moran_seed1.mat';
    randomization_with = 'file';
    HCPtasks_reconstruction_rotateadd;
    save('./results/HCPtasks_reconstruction_200modes_5000rotates_moran.mat','recon_rotateadd','recon_rotate','corr_geom_beta','corr_orig_rotate','zstat_avg');
    recon_acc_n200_rnd5000 = recon_rotateadd;
    recon_acc_n200_rnd5000_no_add = recon_rotate;
    corr_geom_beta_moran = corr_geom_beta;
    corr_orig_rotate_moran = corr_orig_rotate;
    zstat_avg_moran = zstat_avg;
end
