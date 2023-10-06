clear
close all
clc
addpath(genpath('../BrainEigenmodes/functions_matlab'));
addpath(genpath('../rand_index'));

% cortex and medial wall mask
% ---------------------------
cortex = dlmread('../BrainEigenmodes/data/template_surfaces_volumes/fsLR_32k_cortex-lh_mask.txt');
cortex_ind = find(cortex);

% Prepare results in figures
% --------------------------
prepare_results;
recon_acc_without_rotate_n200 = recon_acc_n200([44,19,41,11,47,3,30],:); % the key-tasks

% Read surface for visualization
% ------------------------------
[vertices, faces] = read_vtk(sprintf('../BrainEigenmodes/data/template_surfaces_volumes/fsLR_32k_midthickness-lh.vtk'));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
%% Figure 1

% Parameters
% ----------
nRotation_gambling = 36;     % for the spin test
nRotation_relational = 1288; % for the spin test

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel','Parcel + Connectome'};

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
fig = figure(1);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Figure 1');

% --------------------------
% a: Reconstruction accuracy
% --------------------------
ax1 = axes(fig);
ax1.Units = 'normalized';
ax1.Position = [0.055,0.70,0.3,0.24];
annotation(fig, 'textbox', [0.01, 0.99, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
hold on
for nEig = 1:size(recon_acc_n200,2)
    Y = recon_acc_n200(:,nEig);
    [X] = plot_violin_scatter(Y,nEig,17,0.5,0);
    scatter(ax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(nEig,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(ax1,Y,'Positions',nEig,'Orientation','vertical','Width',0.2,'Colors',colors2(nEig,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(nEig,:)*0.5);
    
    % Label
    text(nEig,0.69,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',45);
end
[p,h,s] = signrank(recon_acc_n200(:,1) - recon_acc_n200(:,6));
if p < 0.05
    line(ax1,[1,1,6,6],[max(recon_acc_n200(:,1))+0.005,1,1,max(recon_acc_n200(:,6))+0.005],'Color','k');
end
X0 = nEig;
for nEig = [1,5,6]
    X0 = X0 + 1;
    Y = recon_acc_n100(:,nEig);
    [X] = plot_violin_scatter(Y,X0,17,0.5,0);
    scatter(ax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(nEig,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(ax1,Y,'Positions',X0,'Orientation','vertical','Width',0.2,'Colors',colors2(nEig,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(nEig,:)*0.5);
    
    % Label
    text(ax1,X0,0.69,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',45);
end
for nEig = [1,5,6]
    X0 = X0 + 1;
    Y = recon_acc_n50(:,nEig);
    [X] = plot_violin_scatter(Y,X0,17,0.5,0);
    scatter(ax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(nEig,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(ax1,Y,'Positions',X0,'Orientation','vertical','Width',0.2,'Colors',colors2(nEig,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(nEig,:)*0.5);
    
    % Label
    text(ax1,X0,0.69,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',45);
end
[p,h,s] = signrank(recon_acc_n100(:,1) - recon_acc_n100(:,6));
[p,h,s] = signrank(recon_acc_n50(:,1) - recon_acc_n50(:,6));
if p < 0.05
    line(ax1,[10,10,12,12],[max(recon_acc_n50(:,1))+0.01,0.965,0.965,max(recon_acc_n50(:,6))+0.01],'Color','k');
end
xlim([0.5,X0+0.5]);
ylim([0.70,1]);
set(ax1,'FontSize',10,'FontWeight','n','XColor','w','box','off');
ax1.Box = 'off';
ylabel('Reconstruction accuracy')
line(ax1,[6.5,6.5],[0.72,0.98],'Color','k','LineStyle',':')
line(ax1,[9.5,9.5],[0.72,0.98],'Color','k','LineStyle',':')
text(ax1,3.5,0.72,'200 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
text(ax1,8,0.72,'100 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
text(ax1,11,0.72,'50 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
ax1.Position = [0.055,0.70,0.3,0.24];

% -------------------------------------
% b: Eigen modes of Schaefer conenctome
% -------------------------------------
cmap = hot(500); cmap = cmap(end:-1:1,:);
annotation(fig, 'textbox', [0.36, 0.99, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
ax_annot1 = axes(fig);
ax_annot1.Position = [0.40,0.89,0.01,0.01];
xlim(ax_annot1,[-1,1]);ylim(ax_annot1,[-1,1]);axis(ax_annot1,'off');
text(ax_annot1,0,0,'Parcel','FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');
ax_annot2 = axes(fig);
ax_annot2.Position = [0.40,0.76,0.01,0.01];
xlim(ax_annot2,[-1,1]);ylim(ax_annot2,[-1,1]);axis(ax_annot2,'off');
text(ax_annot2,0,0,sprintf('Parcel\n+\nConnectome'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

parcels = [23,44,150,163,172,180];
for nMode = 1:6
    pos = [0.45+0.09*(nMode-1),0.83,0.06,0.12];
    eval(sprintf('ax2_%i = axes(fig);',nMode+6));
    eval(sprintf('temp_ax = ax2_%i;',nMode+6));
    temp_ax.Units = 'normalized';
    temp_ax.Position = pos;
    nParcel = parcels(nMode);
    cdata = eigenmodes_n200.Schaefer_mask(:,nParcel); cdata = cdata / max(cdata);
    patch(temp_ax, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
    caxis(temp_ax,[0,1]);
    colormap(temp_ax,cmap);
    material dull
    axis vis3d
    axis off
    axis image
    view(270,0);
    camlight(temp_ax,'headlight')
    title(sprintf('Mode %i',nParcel),'FontSize',10,'FontWeight','n')
    eval(sprintf('ax2_%i = temp_ax;',nMode+6));
end
for nMode = 1:6
    pos = [0.45+0.09*(nMode-1),0.70,0.06,0.12];
    eval(sprintf('ax2_%i = axes(fig);',nMode));
    eval(sprintf('temp_ax = ax2_%i;',nMode));
    temp_ax.Units = 'normalized';
    temp_ax.Position = pos;
    nParcel = parcels(nMode);
    cdata = eigenmodes_n200.Schaefer_Connectome(:,nParcel); cdata = cdata / max(cdata);
    patch(temp_ax, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
    caxis(temp_ax,[0,1]);
    colormap(temp_ax,cmap);
    material dull
    axis vis3d
    axis off
    axis image
    view(270,0);
    camlight(temp_ax,'headlight')
    eval(sprintf('ax2_%i = temp_ax;',nMode+6));
end

% ---------------------
% c: Effect of rotation
% ---------------------
ax3 = axes(fig);
ax3.Units = 'normalized';
ax3.Position = [0.06,0.12,0.19,0.25];
annotation(fig, 'textbox', [0.01, 0.48, 0.01, 0.01], 'string', 'c', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
hold on
line(ax3,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(ax3,[0,0],[-0.4,0.8],'Color','k','LineStyle',':','LineWidth',1);
scatter(ax3,corr_geom_beta_spin.gambling_punish_reward,corr_orig_rotate_spin.gambling_punish_reward,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',0.3)
scatter(ax3,corr_geom_beta_spin.relational_match_rel,corr_orig_rotate_spin.relational_match_rel,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(7,:),'MarkerFaceAlpha',0.3)
xlim(ax3,[-1,1]);
ylim(ax3,[-0.4,0.8]);
xlabel(ax3,'Correlation between \beta values')
ylabel(ax3,'Correlation with rotated maps')

x = corr_geom_beta_spin.gambling_punish_reward;
y = corr_orig_rotate_spin.gambling_punish_reward;
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(ax3,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(3,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(ax3,-0.95,0.74,sprintf('Pearson r = %.03f',r),'Color',colors(3,:),'HorizontalAlignment','left');
r = corr(x',y','type','Spearman');
text(ax3,-0.95,0.64,sprintf('Spearman r = %.03f',r),'Color',colors(3,:),'HorizontalAlignment','left');

x = corr_geom_beta_spin.relational_match_rel;
y = corr_orig_rotate_spin.relational_match_rel;
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(ax3,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(7,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(ax3,-0.95,0.51,sprintf('Pearson r = %.03f',r),'Color',colors(7,:),'HorizontalAlignment','left');
r = corr(x',y','type','Spearman');
text(ax3,-0.95,0.41,sprintf('Spearman r = %.03f',r),'Color',colors(7,:),'HorizontalAlignment','left');

ax3_r = axes(fig);
ax3_r.Units = 'normalized';
ax3_r.Position = [0.255,0.12,0.04,0.25];
Y = corr_orig_rotate_spin.relational_match_rel;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(ax3_r,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_r,Vq,Xq,'-r','Color',colors(7,:));
axis(ax3_r,'off');
ylim(ax3_r,[-0.4,0.8]);

Y = corr_orig_rotate_spin.gambling_punish_reward;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(ax3_r,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_r,Vq,Xq,'-k','Color',colors(3,:))
axis(ax3_r,'off');
ylim(ax3_r,[-0.4,0.8]);

ax3_u = axes(fig);
ax3_u.Units = 'normalized';
ax3_u.Position = [0.06,0.38,0.19,0.08];
X = corr_geom_beta_spin.relational_match_rel;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(ax3_u,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_u,Xq,Vq,'-r','Color',colors(7,:));
axis(ax3_u,'off');
xlim(ax3_u,[-1,1]);

X = corr_geom_beta_spin.gambling_punish_reward;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(ax3_u,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_u,Xq,Vq,'-k','Color',colors(3,:));
axis(ax3_u,'off');
xlim(ax3_u,[-1,1]);

ax3_text1 = axes(fig);
ax3_text1.Units = 'normalized';
ax3_text1.Position = [0.19,0.48,0.01,0.01];
xlim(ax3_text1,[-1,1]);ylim(ax3_text1,[-1,1]);axis(ax3_text1,'off');
text(ax3_text1,0,0,'Gambling','FontSize',10,'FontWeight','n','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','Color',colors(3,:));

ax3_text2 = axes(fig);
ax3_text2.Units = 'normalized';
ax3_text2.Position = [0.19,0.42,0.01,0.01];
xlim(ax3_text2,[-1,1]);ylim(ax3_text2,[-1,1]);axis(ax3_text2,'off');
text(ax3_text2,0,0,'Relational','FontSize',10,'FontWeight','n','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','Color',colors(7,:));

ax3_brain1 = axes(fig);
ax3_brain1.Units = 'normalized';
ax3_brain1.Position = [0.25,0.45,0.03,0.06];
cdata = zstat_avg_spin.gambling_punish_reward;
patch(ax3_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain1,bluewhitered);
material(ax3_brain1,'dull');
axis(ax3_brain1,'vis3d');
axis(ax3_brain1,'off');
axis(ax3_brain1,'image');
view(ax3_brain1,270,0);
camlight(ax3_brain1,'headlight')
title(ax3_brain1,'Original','FontSize',10,'FontWeight','n');

ax3_brain2 = axes(fig);
ax3_brain2.Units = 'normalized';
ax3_brain2.Position = [0.30,0.45,0.03,0.06];
cdata = squeeze(zstat_rotates.gambling_punish_reward{1,1}(:,1,nRotation_gambling));
patch(ax3_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain2,bluewhitered);
material(ax3_brain2,'dull')
axis(ax3_brain2,'vis3d')
axis(ax3_brain2,'off');
axis(ax3_brain2,'image');
view(ax3_brain2,270,0);
camlight(ax3_brain2,'headlight')
scatter(ax3,corr_geom_beta_spin.gambling_punish_reward(1,nRotation_gambling),corr_orig_rotate_spin.gambling_punish_reward(1,nRotation_gambling),20,'Marker','o','MarkerEdgeColor',colors(3,:)*0.6,'LineWidth',2);
ax3.Position = [0.06,0.12,0.20,0.26];
title(ax3_brain2,'Rotated','FontSize',10,'FontWeight','n');

ax3_brain3 = axes(fig);
ax3_brain3.Units = 'normalized';
ax3_brain3.Position = [0.25,0.39,0.03,0.06];
cdata = zstat_avg_spin.relational_match_rel;
patch(ax3_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain3,bluewhitered);
material(ax3_brain3,'dull');
axis(ax3_brain3,'vis3d');
axis(ax3_brain3,'off');
axis(ax3_brain3,'image');
view(ax3_brain3,270,0);
camlight(ax3_brain3,'headlight')

ax3_brain4 = axes(fig);
ax3_brain4.Units = 'normalized';
ax3_brain4.Position = [0.30,0.39,0.03,0.06];
cdata = squeeze(zstat_rotates.relational_match_rel{1,1}(:,1,nRotation_relational));
patch(ax3_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain4,bluewhitered);
material(ax3_brain4,'dull')
axis(ax3_brain4,'vis3d')
axis(ax3_brain4,'off');
axis(ax3_brain4,'image');
view(ax3_brain4,270,0);
camlight(ax3_brain4,'headlight')
scatter(ax3,corr_geom_beta_spin.relational_match_rel(1,nRotation_relational),corr_orig_rotate_spin.relational_match_rel(1,nRotation_relational),20,'Marker','o','MarkerEdgeColor',colors(7,:)*0.6,'LineWidth',2);

ax3_accu = axes(fig);
ax3_accu.Units = 'normalized';
ax3_accu.Position = [0.30,0.26,0.03,0.12];
hold on
Y1 = recon_acc_n200_rot5000_no_add.gambling_punish_reward(:,1);
Y2 = recon_acc_n200_rot5000_no_add.relational_match_rel(:,1);
val_rot_no_add = Y1;
[X] = plot_violin_scatter(val_rot_no_add,1,40,0.5,0);
scatter(ax3_accu,X,val_rot_no_add,5,'Marker','o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
hbp1 = boxplot(ax3_accu,val_rot_no_add,'Positions',1,'Orientation','vertical','Width',0.3,'Colors',colors(3,:)*0.5,'outliersize',1.5,'symbol','k+');
set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors(3,:)*0.5);
val_rot_no_add = Y2;
[X] = plot_violin_scatter(val_rot_no_add,2,40,0.5,0);
scatter(ax3_accu,X,val_rot_no_add,5,'Marker','o','MarkerFaceColor',colors(7,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
hbp1 = boxplot(ax3_accu,val_rot_no_add,'Positions',2,'Orientation','vertical','Width',0.3,'Colors',colors(7,:)*0.5,'outliersize',1.5,'symbol','k+');
set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors(7,:)*0.5);
xlim([0.5,2.5]);
ylim([0.9,1]);
axis(ax3_accu,'on');
set(ax3_accu,'FontSize',10,'FontWeight','n','XColor','w','box','off','YTick',[0.9,1.0]);
ylabel(ax3_accu,sprintf('Recon. acc.'));
text(ax3_accu,1.5,0.88,'Geometric','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
line(ax3_accu,[1,1,2,2],[0.994,1,1,0.997],'Color','k');
ax3_accu.Position = [0.30,0.26,0.03,0.12];

ax3.Position = [0.06,0.12,0.19,0.25];

% ---------------------------------------------------------------------------
% d: Reconstruction accuracy with rotated task activity maps (5000 rotations)
% ---------------------------------------------------------------------------
% Task labels
% -----------
str_tasks = {'emotion_faces_shapes','gambling_punish_reward','language_math_story','motor_cue_avg','relational_match_rel','social_tom_random','wm_2bk_0bk'};
str_tasks_fig = {'Emotion','Gambling','Language','Motor','Relational','Social',sprintf('Working\nmemory')};
ax4 = axes(fig);
ax4.Units = 'normalized';
ax4.Position = [0.40,0.13,0.60,0.55];
annotation(fig, 'textbox', [0.36, 0.71, 0.01, 0.01], 'string', 'd', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
hold on
for nTask = 1:7
    X0 = 0;
    for nMode = [1,2,5,6]
        X0 = X0 + 1;
        X1 = nTask-0.24+(X0-1)*0.16;
        val_no_rot = recon_acc_without_rotate_n200(nTask,nMode);
        eval(sprintf('val_rot    = recon_acc_n200_rot5000_no_add.%s(:,nMode);',str_tasks{nTask}));
        [X] = plot_violin_scatter(val_rot,X1,40,0.1,0);
        scatter(ax4,X,val_rot,5,'Marker','o','MarkerFaceColor',colors2(nMode,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
        hbp1 = boxplot(ax4,val_rot,'Positions',X1,'Orientation','vertical','Width',0.07,'Colors',colors2(nMode,:)*0.5,'outliersize',1.5,'symbol','k+');
        set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors2(nMode,:)*0.5);
        line([X1-0.1,X1+0.1],[1,1]*val_no_rot,'Color',colors2(3,:)*0.8,'LineWidth',1,'LineStyle','-');
        [p,h,s] = signrank(val_rot - val_no_rot);
        p = 1 - (1 - p)^7; % Bonferroni correction
        e = s.zval/sqrt(length(val_rot));
        m = max(val_rot);
        if p < 0.05 && p >=0.0001
            if median(val_rot) < val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','none');
            elseif median(val_rot) > val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','none');
            end
        elseif p < 0.0001
            if median(val_rot) < val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k');
            elseif median(val_rot) > val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k');
            end
        end
        disp(e)
    end
    
    % Task name
    text(nTask,0.79,str_tasks_fig{nTask},'FontSize',10,'FontWeight','n','HorizontalAlignment','center');
end
xlim([0.5,7.5]);
ylim([0.80,1]);
set(ax4,'FontSize',10,'FontWeight','n','XColor','w','box','off');
ax4.Box = 'off';
ylabel('Reconstruction accuracy')
scatter(ax4,3.15,0.885,15,'Marker','o','MarkerFaceColor',colors2(1,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.875,15,'Marker','o','MarkerFaceColor',colors2(2,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.865,15,'Marker','o','MarkerFaceColor',colors2(5,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.855,15,'Marker','o','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
line(ax4,[3.12,3.18],[1,1]*0.845,'Color',colors2(3,:)*0.8,'LineWidth',1,'LineStyle','-');
scatter(ax4,3.15,0.835,15,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(ax4,3.15,0.825,15,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(ax4,3.15,0.815,12,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(ax4,3.15,0.805,12,'Marker','+','MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',1.5);
text(ax4,3.2,0.885,'Geometric','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.875,'Connectome','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.865,'Parcel','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.855,'Parcel + Connectome','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.845,'Original','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.835,'Corrected \itp\rm < 0.0001 (overfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.825,'Corrected \itp\rm < 0.0001 (underfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.815,'Corrected \itp\rm < 0.05 (underfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.805,'Outliers','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
ax4.Position = [0.40,0.13,0.60,0.55];

ax4_brain1 = axes(fig);
ax4_brain1.Units = 'normalized';
ax4_brain1.Position = [0.427,0.02,0.03,0.06];
cdata = zstat_avg_spin.emotion_faces_shapes;
patch(ax4_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain1,bluewhitered);
material(ax4_brain1,'dull');
axis(ax4_brain1,'vis3d');
axis(ax4_brain1,'off');
axis(ax4_brain1,'image');
view(ax4_brain1,270,0);
camlight(ax4_brain1,'headlight')

ax4_brain2 = axes(fig);
ax4_brain2.Units = 'normalized';
ax4_brain2.Position = [0.513,0.02,0.03,0.06];
cdata = zstat_avg_spin.gambling_punish_reward;
patch(ax4_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain2,bluewhitered);
material(ax4_brain2,'dull');
axis(ax4_brain2,'vis3d');
axis(ax4_brain2,'off');
axis(ax4_brain2,'image');
view(ax4_brain2,270,0);
camlight(ax4_brain2,'headlight')

ax4_brain3 = axes(fig);
ax4_brain3.Units = 'normalized';
ax4_brain3.Position = [0.599,0.02,0.03,0.06];
cdata = zstat_avg_spin.language_math_story;
patch(ax4_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain3,bluewhitered);
material(ax4_brain3,'dull');
axis(ax4_brain3,'vis3d');
axis(ax4_brain3,'off');
axis(ax4_brain3,'image');
view(ax4_brain3,270,0);
camlight(ax4_brain3,'headlight')

ax4_brain4 = axes(fig);
ax4_brain4.Units = 'normalized';
ax4_brain4.Position = [0.685,0.02,0.03,0.06];
cdata = zstat_avg_spin.motor_cue_avg;
patch(ax4_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain4,bluewhitered);
material(ax4_brain4,'dull');
axis(ax4_brain4,'vis3d');
axis(ax4_brain4,'off');
axis(ax4_brain4,'image');
view(ax4_brain4,270,0);
camlight(ax4_brain4,'headlight')

ax4_brain5 = axes(fig);
ax4_brain5.Units = 'normalized';
ax4_brain5.Position = [0.771,0.02,0.03,0.06];
cdata = zstat_avg_spin.relational_match_rel;
patch(ax4_brain5, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain5,bluewhitered);
material(ax4_brain5,'dull');
axis(ax4_brain5,'vis3d');
axis(ax4_brain5,'off');
axis(ax4_brain5,'image');
view(ax4_brain5,270,0);
camlight(ax4_brain5,'headlight')

ax4_brain6 = axes(fig);
ax4_brain6.Units = 'normalized';
ax4_brain6.Position = [0.857,0.02,0.03,0.06];
cdata = zstat_avg_spin.social_tom_random;
patch(ax4_brain6, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain6,bluewhitered);
material(ax4_brain6,'dull');
axis(ax4_brain6,'vis3d');
axis(ax4_brain6,'off');
axis(ax4_brain6,'image');
view(ax4_brain6,270,0);
camlight(ax4_brain6,'headlight')

ax4_brain7 = axes(fig);
ax4_brain7.Units = 'normalized';
ax4_brain7.Position = [0.943,0.02,0.03,0.06];
cdata = zstat_avg_spin.wm_2bk_0bk;
patch(ax4_brain7, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain7,bluewhitered);
material(ax4_brain7,'dull');
axis(ax4_brain7,'vis3d');
axis(ax4_brain7,'off');
axis(ax4_brain7,'image');
view(ax4_brain7,270,0);
camlight(ax4_brain7,'headlight')

% Save
set(fig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',fig.Position([3,4]) + 1);drawnow;
% saveas(fig,'./figure1.pdf');
%% Supplementary figure 1

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel','Parcel + Connectome'};

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
sfig = figure(2);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 1');

% ---------------------------------
% Figure 1: Reconstruction accuracy
% ---------------------------------
sfax1 = axes(sfig);
sfax1.Units = 'normalized';
sfax1.Position = [0.055,0.67,0.22,0.26];
hold on
for nEig = 1:size(recon_acc_n100,2)
    Y = recon_acc_n100(:,nEig);
    [X] = plot_violin_scatter(Y,nEig,17,0.5,0);
    scatter(sfax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(nEig,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(sfax1,Y,'Positions',nEig,'Orientation','vertical','Width',0.2,'Colors',colors2(nEig,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(nEig,:)*0.5);
    
    % Label
    text(nEig,0.39,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',55);
end
X0 = nEig;
for nEig = 1:size(recon_acc_n50,2)
    X0 = X0 + 1;
    Y = recon_acc_n50(:,nEig);
    [X] = plot_violin_scatter(Y,X0,17,0.5,0);
    scatter(sfax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(nEig,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(sfax1,Y,'Positions',X0,'Orientation','vertical','Width',0.2,'Colors',colors2(nEig,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(nEig,:)*0.5);
    
    % Label
    text(sfax1,X0,0.39,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',55);
end
xlim(sfax1,[0.5,X0+0.5]);
ylim(sfax1,[0.40,1]);
set(sfax1,'FontSize',10,'FontWeight','n','XColor','w','box','off');
sfax1.Box = 'off';
ylabel('Reconstruction accuracy')
line(sfax1,[6.5,6.5],[0.42,0.98],'Color','k','LineStyle',':')
text(sfax1,3.5,0.45,'100 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
text(sfax1,9.5,0.45,'50 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
sfax1.Position = [0.055,0.67,0.22,0.26];

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure1.pdf');
%% Supplementary figure 2

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel (Schaefer)','Parcel (kmedoids)','Parcel (Schaefer) + Connectome','Parcel (kmedoids) + Connectome'};

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
sfig = figure(3);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 2');

% --------------------------
% a: Reconstruction accuracy
% --------------------------
sfax1 = axes(sfig);
sfax1.Units = 'normalized';
sfax1.Position = [0.06,0.67,0.55,0.26];
annotation(sfig, 'textbox', [0.01, 0.99, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
hold(sfax1,'on');
for nEig = 1:8
    if nEig == 5
        Y = recon_acc_n200(:,5);
        color_ind = 5;
    elseif nEig == 6
        Y = recon_acc_kmedoids_n200(:,5);
        color_ind = 5;
    elseif nEig == 7
        Y = recon_acc_n200(:,6);
        color_ind = 6;
    elseif nEig == 8
        Y = recon_acc_kmedoids_n200(:,6);
        color_ind = 6;
    else
        Y = recon_acc_n200(:,nEig);
        color_ind = nEig;
    end
    [X] = plot_violin_scatter(Y,nEig,17,0.5,0);
    scatter(sfax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(color_ind,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(sfax1,Y,'Positions',nEig,'Orientation','vertical','Width',0.2,'Colors',colors2(color_ind,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(color_ind,:)*0.5);
    
    % Label
    text(sfax1,nEig,0.39,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',60);
end
X0 = nEig;
for nEig = 1:8
    X0 = X0 + 1;
    if nEig == 5
        Y = recon_acc_n100(:,5);
        color_ind = 5;
    elseif nEig == 6
        Y = recon_acc_kmedoids_n100(:,5);
        color_ind = 5;
    elseif nEig == 7
        Y = recon_acc_n100(:,6);
        color_ind = 6;
    elseif nEig == 8
        Y = recon_acc_kmedoids_n100(:,6);
        color_ind = 6;
    else
        Y = recon_acc_n100(:,nEig);
        color_ind = nEig;
    end
    [X] = plot_violin_scatter(Y,X0,17,0.5,0);
    scatter(sfax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(color_ind,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(sfax1,Y,'Positions',X0,'Orientation','vertical','Width',0.2,'Colors',colors2(color_ind,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(color_ind,:)*0.5);
    
    % Label
    text(sfax1,X0,0.39,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',60);
end
for nEig = 1:8
    X0 = X0 + 1;
    if nEig == 5
        Y = recon_acc_n50(:,5);
        color_ind = 5;
    elseif nEig == 6
        Y = recon_acc_kmedoids_n50(:,5);
        color_ind = 5;
    elseif nEig == 7
        Y = recon_acc_n50(:,6);
        color_ind = 6;
    elseif nEig == 8
        Y = recon_acc_kmedoids_n50(:,6);
        color_ind = 6;
    else
        Y = recon_acc_n50(:,nEig);
        color_ind = nEig;
    end
    [X] = plot_violin_scatter(Y,X0,17,0.5,0);
    scatter(sfax1,X,Y,5,'Marker','o','MarkerFaceColor',colors2(color_ind,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(sfax1,Y,'Positions',X0,'Orientation','vertical','Width',0.2,'Colors',colors2(color_ind,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(nEig,:)*0.5);
    
    % Label
    text(sfax1,X0,0.39,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',60);
end
xlim(sfax1,[0.5,X0+0.5]);
ylim(sfax1,[0.40,1]);
set(sfax1,'FontSize',10,'FontWeight','n','XColor','w','box','off');
sfax1.Box = 'off';
ylabel('Reconstruction accuracy')
line(sfax1,[8.5,8.5],[0.42,0.98],'Color','k','LineStyle',':')
line(sfax1,[16.5,16.5],[0.42,0.98],'Color','k','LineStyle',':')
text(sfax1,4.5,0.45,'200 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
text(sfax1,12.5,0.45,'100 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
text(sfax1,20.5,0.45,'50 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
sfax1.Position = [0.06,0.67,0.55,0.26];

% -------------
% b: atlas maps
% -------------
annotation(sfig, 'textbox', [0.01, 0.40, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

% surface projection
% ------------------
ax2_brain1 = axes(sfig);
ax2_brain1.Units = 'normalized';
ax2_brain1.Position = [0.08,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(200)];
cdata = schaefer400; cdata(isnan(cdata)) = 0;
patch(ax2_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain1,[0,200]);
colormap(ax2_brain1,cmap);
material(ax2_brain1,'dull');
axis(ax2_brain1,'vis3d');
axis(ax2_brain1,'off');
axis(ax2_brain1,'image');
view(ax2_brain1,270,0);
camlight(ax2_brain1,'headlight')

ax2_brain2 = axes(sfig);
ax2_brain2.Units = 'normalized';
ax2_brain2.Position = [0.17,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(200)];
cdata = schaefer400; cdata(isnan(cdata)) = 0;
patch(ax2_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain2,[0,200]);
colormap(ax2_brain2,cmap);
material(ax2_brain2,'dull');
axis(ax2_brain2,'vis3d');
axis(ax2_brain2,'off');
axis(ax2_brain2,'image');
view(ax2_brain2,90,0);
camlight(ax2_brain2,'headlight')

ax2_brain3 = axes(sfig);
ax2_brain3.Units = 'normalized';
ax2_brain3.Position = [0.26,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(100)];
cdata = schaefer200; cdata(isnan(cdata)) = 0;
patch(ax2_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain3,[0,100]);
colormap(ax2_brain3,cmap);
material(ax2_brain3,'dull');
axis(ax2_brain3,'vis3d');
axis(ax2_brain3,'off');
axis(ax2_brain3,'image');
view(ax2_brain3,270,0);
camlight(ax2_brain3,'headlight')

ax2_brain4 = axes(sfig);
ax2_brain4.Units = 'normalized';
ax2_brain4.Position = [0.35,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(100)];
cdata = schaefer200; cdata(isnan(cdata)) = 0;
patch(ax2_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain4,[0,100]);
colormap(ax2_brain4,cmap);
material(ax2_brain4,'dull');
axis(ax2_brain4,'vis3d');
axis(ax2_brain4,'off');
axis(ax2_brain4,'image');
view(ax2_brain4,90,0);
camlight(ax2_brain4,'headlight')

ax2_brain5 = axes(sfig);
ax2_brain5.Units = 'normalized';
ax2_brain5.Position = [0.44,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(50)];
cdata = schaefer100; cdata(isnan(cdata)) = 0;
patch(ax2_brain5, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain5,[0,100]);
colormap(ax2_brain5,cmap);
material(ax2_brain5,'dull');
axis(ax2_brain5,'vis3d');
axis(ax2_brain5,'off');
axis(ax2_brain5,'image');
view(ax2_brain5,270,0);
camlight(ax2_brain5,'headlight')

ax2_brain6 = axes(sfig);
ax2_brain6.Units = 'normalized';
ax2_brain6.Position = [0.53,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(50)];
cdata = schaefer100; cdata(isnan(cdata)) = 0;
patch(ax2_brain6, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain6,[0,50]);
colormap(ax2_brain6,cmap);
material(ax2_brain6,'dull');
axis(ax2_brain6,'vis3d');
axis(ax2_brain6,'off');
axis(ax2_brain6,'image');
view(ax2_brain6,90,0);
camlight(ax2_brain6,'headlight')

ax_annot1 = axes(sfig);
ax_annot1.Position = [0.03,0.25,0.01,0.01];
xlim(ax_annot1,[-1,1]);ylim(ax_annot1,[-1,1]);axis(ax_annot1,'off');
text(ax_annot1,0,0,sprintf('Schaefer'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax_annot2 = axes(sfig);
ax_annot2.Position = [0.15,0.34,0.01,0.01];
xlim(ax_annot2,[-1,1]);ylim(ax_annot2,[-1,1]);axis(ax_annot2,'off');
text(ax_annot2,0,0,sprintf('200 parcels'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax_annot3 = axes(sfig);
ax_annot3.Position = [0.33,0.34,0.01,0.01];
xlim(ax_annot3,[-1,1]);ylim(ax_annot3,[-1,1]);axis(ax_annot3,'off');
text(ax_annot3,0,0,sprintf('100 parcels'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax_annot4 = axes(sfig);
ax_annot4.Position = [0.51,0.34,0.01,0.01];
xlim(ax_annot4,[-1,1]);ylim(ax_annot4,[-1,1]);axis(ax_annot4,'off');
text(ax_annot4,0,0,sprintf('50 parcels'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax_annot5 = axes(sfig);
ax_annot5.Position = [0.03,0.10,0.01,0.01];
xlim(ax_annot5,[-1,1]);ylim(ax_annot5,[-1,1]);axis(ax_annot5,'off');
text(ax_annot5,0,0,sprintf('kmedoids'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax2_brain7 = axes(sfig);
ax2_brain7.Units = 'normalized';
ax2_brain7.Position = [0.08,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(200)];
cdata = atlasfree200; cdata(isnan(cdata)) = 0;
patch(ax2_brain7, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain7,[0,200]);
colormap(ax2_brain7,cmap);
material(ax2_brain7,'dull');
axis(ax2_brain7,'vis3d');
axis(ax2_brain7,'off');
axis(ax2_brain7,'image');
view(ax2_brain7,270,0);
camlight(ax2_brain7,'headlight')

ax2_brain8 = axes(sfig);
ax2_brain8.Units = 'normalized';
ax2_brain8.Position = [0.17,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(200)];
cdata = atlasfree200; cdata(isnan(cdata)) = 0;
patch(ax2_brain8, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain8,[0,200]);
colormap(ax2_brain8,cmap);
material(ax2_brain8,'dull');
axis(ax2_brain8,'vis3d');
axis(ax2_brain8,'off');
axis(ax2_brain8,'image');
view(ax2_brain8,90,0);
camlight(ax2_brain8,'headlight')

ax2_brain9 = axes(sfig);
ax2_brain9.Units = 'normalized';
ax2_brain9.Position = [0.26,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(100)];
cdata = atlasfree100; cdata(isnan(cdata)) = 0;
patch(ax2_brain9, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain9,[0,100]);
colormap(ax2_brain9,cmap);
material(ax2_brain9,'dull');
axis(ax2_brain9,'vis3d');
axis(ax2_brain9,'off');
axis(ax2_brain9,'image');
view(ax2_brain9,270,0);
camlight(ax2_brain9,'headlight')

ax2_brain10 = axes(sfig);
ax2_brain10.Units = 'normalized';
ax2_brain10.Position = [0.35,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(100)];
cdata = atlasfree100; cdata(isnan(cdata)) = 0;
patch(ax2_brain10, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain10,[0,100]);
colormap(ax2_brain10,cmap);
material(ax2_brain10,'dull');
axis(ax2_brain10,'vis3d');
axis(ax2_brain10,'off');
axis(ax2_brain10,'image');
view(ax2_brain10,90,0);
camlight(ax2_brain10,'headlight')

ax2_brain11 = axes(sfig);
ax2_brain11.Units = 'normalized';
ax2_brain11.Position = [0.44,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(50)];
cdata = atlasfree50; cdata(isnan(cdata)) = 0;
patch(ax2_brain11, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain11,[0,100]);
colormap(ax2_brain11,cmap);
material(ax2_brain11,'dull');
axis(ax2_brain11,'vis3d');
axis(ax2_brain11,'off');
axis(ax2_brain11,'image');
view(ax2_brain11,270,0);
camlight(ax2_brain11,'headlight')

ax2_brain12 = axes(sfig);
ax2_brain12.Units = 'normalized';
ax2_brain12.Position = [0.53,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(50)];
cdata = atlasfree50; cdata(isnan(cdata)) = 0;
patch(ax2_brain12, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax2_brain12,[0,50]);
colormap(ax2_brain12,cmap);
material(ax2_brain12,'dull');
axis(ax2_brain12,'vis3d');
axis(ax2_brain12,'off');
axis(ax2_brain12,'image');
view(ax2_brain12,90,0);
camlight(ax2_brain12,'headlight')

ax_annot6 = axes(sfig);
ax_annot6.Position = [0.16,0.18,0.01,0.01];
xlim(ax_annot6,[-1,1]);ylim(ax_annot6,[-1,1]);axis(ax_annot6,'off');
text(ax_annot6,0,0,sprintf('Adjusted Rand index = %.03f',ari200),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','right','VerticalAlignment','middle');

ax_annot7 = axes(sfig);
ax_annot7.Position = [0.34,0.18,0.01,0.01];
xlim(ax_annot7,[-1,1]);ylim(ax_annot7,[-1,1]);axis(ax_annot7,'off');
text(ax_annot7,0,0,sprintf('%.03f',ari100),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','right','VerticalAlignment','middle');

ax_annot8 = axes(sfig);
ax_annot8.Position = [0.52,0.18,0.01,0.01];
xlim(ax_annot8,[-1,1]);ylim(ax_annot8,[-1,1]);axis(ax_annot8,'off');
text(ax_annot8,0,0,sprintf('%.03f',ari50),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','right','VerticalAlignment','middle');

fprintf('\nGlasser vs. Schaefer200 = %.03f\n',ari_gla_sch200)
fprintf('Glasser vs. Schaefer100 = %.03f\n',ari_gla_sch100)
fprintf('Glasser vs. Schaefer50  = %.03f\n\n',ari_gla_sch50)
fprintf('Glasser vs. kmedoids200 = %.03f\n',ari_gla_kmd200)
fprintf('Glasser vs. kmedoids100 = %.03f\n',ari_gla_kmd100)
fprintf('Glasser vs. kmedoids50  = %.03f\n',ari_gla_kmd50)

% ------------------------------------
% c: Reconstruction accuracy (Glasser)
% ------------------------------------
sfax3 = axes(sfig);
sfax3.Units = 'normalized';
sfax3.Position = [0.66,0.67,0.184,0.26];
annotation(sfig, 'textbox', [0.63, 0.99, 0.01, 0.01], 'string', 'c', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
hold(sfax3,'on');
str_modes{5} = 'Parcel (Glasser)';
str_modes{7} = 'Parcel (Glasser) + Connectome';
for nEig = 1:8
    if nEig == 5
        Y = recon_acc_n180(:,5);
        color_ind = 5;
    elseif nEig == 6
        Y = recon_acc_kmedoids_n180(:,5);
        color_ind = 5;
    elseif nEig == 7
        Y = recon_acc_n180(:,6);
        color_ind = 6;
    elseif nEig == 8
        Y = recon_acc_kmedoids_n180(:,6);
        color_ind = 6;
    else
        Y = recon_acc_n180(:,nEig);
        color_ind = nEig;
    end
    [X] = plot_violin_scatter(Y,nEig,17,0.5,0);
    scatter(sfax3,X,Y,5,'Marker','o','MarkerFaceColor',colors2(color_ind,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    hbp = boxplot(sfax3,Y,'Positions',nEig,'Orientation','vertical','Width',0.2,'Colors',colors2(color_ind,:)*0.5,'outliersize',1.5,'symbol','k+');
    set(hbp,{'linew'},{1});set(hbp(6),'Color',colors2(color_ind,:)*0.5);
    
    % Label
    text(sfax3,nEig,0.39,str_modes{nEig},'FontSize',10,'FontWeight','n','HorizontalAlignment','right','Rotation',60);
end
X0 = nEig;
xlim(sfax3,[0.5,X0+0.5]);
ylim(sfax3,[0.40,1]);
set(sfax3,'FontSize',10,'FontWeight','n','XColor','w','box','off');
sfax3.Box = 'off';
text(sfax3,4.5,0.45,'180 modes','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
sfax3.Position = [0.66,0.67,0.184,0.26];

% -------------
% d: atlas maps
% -------------
annotation(sfig, 'textbox', [0.63, 0.40, 0.01, 0.01], 'string', 'd', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

ax4_brain1 = axes(sfig);
ax4_brain1.Units = 'normalized';
ax4_brain1.Position = [0.69,0.05,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(180)];
cdata = atlasfree180; cdata(isnan(cdata)) = 0;
patch(ax4_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax4_brain1,[0,180]);
colormap(ax4_brain1,cmap);
material(ax4_brain1,'dull');
axis(ax4_brain1,'vis3d');
axis(ax4_brain1,'off');
axis(ax4_brain1,'image');
view(ax4_brain1,270,0);
camlight(ax4_brain1,'headlight')

ax4_brain2 = axes(sfig);
ax4_brain2.Units = 'normalized';
ax4_brain2.Position = [0.78,0.05,0.06,0.12];
patch(ax4_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax4_brain2,[0,180]);
colormap(ax4_brain2,cmap);
material(ax4_brain2,'dull');
axis(ax4_brain2,'vis3d');
axis(ax4_brain2,'off');
axis(ax4_brain2,'image');
view(ax4_brain2,90,0);
camlight(ax4_brain2,'headlight')

ax4_brain3 = axes(sfig);
ax4_brain3.Units = 'normalized';
ax4_brain3.Position = [0.69,0.20,0.06,0.12];
cmap = [0.5,0.5,0.5;lines(180)];
cdata = glasser360-180; cdata(isnan(cdata)) = 0;
patch(ax4_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax4_brain3,[0,180]);
colormap(ax4_brain3,cmap);
material(ax4_brain3,'dull');
axis(ax4_brain3,'vis3d');
axis(ax4_brain3,'off');
axis(ax4_brain3,'image');
view(ax4_brain3,270,0);
camlight(ax4_brain3,'headlight')

ax4_brain4 = axes(sfig);
ax4_brain4.Units = 'normalized';
ax4_brain4.Position = [0.78,0.20,0.06,0.12];
patch(ax4_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis(ax4_brain4,[0,180]);
colormap(ax4_brain4,cmap);
material(ax4_brain4,'dull');
axis(ax4_brain4,'vis3d');
axis(ax4_brain4,'off');
axis(ax4_brain4,'image');
view(ax4_brain4,90,0);
camlight(ax4_brain4,'headlight')

ax_annot9 = axes(sfig);
ax_annot9.Position = [0.77,0.18,0.01,0.01];
xlim(ax_annot9,[-1,1]);ylim(ax_annot9,[-1,1]);axis(ax_annot9,'off');
text(ax_annot9,0,0,sprintf('%.03f',ari180),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','right','VerticalAlignment','middle');

ax_annot10 = axes(sfig);
ax_annot10.Position = [0.64,0.25,0.01,0.01];
xlim(ax_annot10,[-1,1]);ylim(ax_annot10,[-1,1]);axis(ax_annot10,'off');
text(ax_annot10,0,0,sprintf('Glasser'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax_annot11 = axes(sfig);
ax_annot11.Position = [0.76,0.34,0.01,0.01];
xlim(ax_annot11,[-1,1]);ylim(ax_annot11,[-1,1]);axis(ax_annot11,'off');
text(ax_annot11,0,0,sprintf('180 parcels'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

ax_annot12 = axes(sfig);
ax_annot12.Position = [0.64,0.10,0.01,0.01];
xlim(ax_annot12,[-1,1]);ylim(ax_annot12,[-1,1]);axis(ax_annot12,'off');
text(ax_annot12,0,0,sprintf('kmedoids'),'FontSize',10,'FontWeight','n','Rotation',00,'HorizontalAlignment','center','VerticalAlignment','middle');

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure2.pdf');
%% Supplementary figure 3

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel','Parcel + Connectome'};

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
sfig = figure(4);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 3');

% Social
% ------
sfax2_1 = axes(sfig);
sfax2_1.Units = 'normalized';
sfax2_1.Position = [0.34,0.75,0.10,0.15];
sfax2_1r = axes(sfig);
sfax2_1r.Units = 'normalized';
sfax2_1r.Position = [0.445,0.75,0.02,0.15];
sfax2_1u = axes(sfig);
sfax2_1u.Units = 'normalized';
sfax2_1u.Position = [0.34,0.91,0.10,0.04];
hold(sfax2_1,'on');
hold(sfax2_1r,'on');
hold(sfax2_1u,'on');

line(sfax2_1,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_1,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_1,corr_geom_beta_spin.social_tom_random,corr_orig_rotate_spin.social_tom_random,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(1,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_1,[-1,1]); ylim(sfax2_1,[-1,1]);
ylabel(sfax2_1,sprintf('Correlation with\nrotated maps'))
x = corr_geom_beta_spin.social_tom_random;
y = corr_orig_rotate_spin.social_tom_random;
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_1,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(1,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_1,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(1,:),'HorizontalAlignment','left');

Y = corr_orig_rotate_spin.social_tom_random;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_1r,BINS,N,'FaceColor',colors(1,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_1r,Vq,Xq,'-r','Color',colors(1,:));
axis(sfax2_1r,'off');
ylim(sfax2_1r,[-1,1]);

X = corr_geom_beta_spin.social_tom_random;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_1u,BINS,N,'FaceColor',colors(1,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_1u,Xq,Vq,'-r','Color',colors(1,:));
axis(sfax2_1u,'off');
xlim(sfax2_1u,[-1,1]);
title(sfax2_1u,'Social','FontSize',10,'FontWeight','n');

% Motor
% -----
sfax2_2 = axes(sfig);
sfax2_2.Units = 'normalized';
sfax2_2.Position = [0.50,0.75,0.10,0.15];
sfax2_2r = axes(sfig);
sfax2_2r.Units = 'normalized';
sfax2_2r.Position = [0.605,0.75,0.02,0.15];
sfax2_2u = axes(sfig);
sfax2_2u.Units = 'normalized';
sfax2_2u.Position = [0.50,0.91,0.10,0.04];
hold(sfax2_2,'on');
hold(sfax2_2r,'on');
hold(sfax2_2u,'on');

line(sfax2_2,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_2,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
x = corr_geom_beta_spin.motor_cue_avg;
y = corr_orig_rotate_spin.motor_cue_avg;
scatter(sfax2_2,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_2,[-1,1]); ylim(sfax2_2,[-1,1]);
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_2,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(2,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_2,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(2,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_2r,BINS,N,'FaceColor',colors(2,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_2r,Vq,Xq,'-r','Color',colors(2,:));
axis(sfax2_2r,'off');
ylim(sfax2_2r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_2u,BINS,N,'FaceColor',colors(2,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_2u,Xq,Vq,'-r','Color',colors(2,:));
axis(sfax2_2u,'off');
xlim(sfax2_2u,[-1,1]);
title(sfax2_2u,'Motor','FontSize',10,'FontWeight','n');

% Gambling
% --------
sfax2_3 = axes(sfig);
sfax2_3.Units = 'normalized';
sfax2_3.Position = [0.66,0.75,0.10,0.15];
sfax2_3r = axes(sfig);
sfax2_3r.Units = 'normalized';
sfax2_3r.Position = [0.765,0.75,0.02,0.15];
sfax2_3u = axes(sfig);
sfax2_3u.Units = 'normalized';
sfax2_3u.Position = [0.66,0.91,0.10,0.04];
hold(sfax2_3,'on');
hold(sfax2_3r,'on');
hold(sfax2_3u,'on');

line(sfax2_3,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_3,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
x = corr_geom_beta_spin.gambling_punish_reward;
y = corr_orig_rotate_spin.gambling_punish_reward;
scatter(sfax2_3,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_3,[-1,1]); ylim(sfax2_3,[-1,1]);
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_3,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(3,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_3,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(3,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_3r,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_3r,Vq,Xq,'-r','Color',colors(3,:));
axis(sfax2_3r,'off');
ylim(sfax2_3r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_3u,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_3u,Xq,Vq,'-r','Color',colors(3,:));
axis(sfax2_3u,'off');
xlim(sfax2_3u,[-1,1]);
title(sfax2_3u,'Gambling','FontSize',10,'FontWeight','n');

% Working memory
% --------------
sfax2_4 = axes(sfig);
sfax2_4.Units = 'normalized';
sfax2_4.Position = [0.82,0.75,0.10,0.15];
sfax2_4r = axes(sfig);
sfax2_4r.Units = 'normalized';
sfax2_4r.Position = [0.925,0.75,0.02,0.15];
sfax2_4u = axes(sfig);
sfax2_4u.Units = 'normalized';
sfax2_4u.Position = [0.82,0.91,0.10,0.04];
hold(sfax2_4,'on');
hold(sfax2_4r,'on');
hold(sfax2_4u,'on');

line(sfax2_4,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_4,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
x = corr_geom_beta_spin.wm_2bk_0bk;
y = corr_orig_rotate_spin.wm_2bk_0bk;
scatter(sfax2_4,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(4,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_4,[-1,1]); ylim(sfax2_4,[-1,1]);
xlabel(sfax2_4,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_4,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(4,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_4,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(4,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_4r,BINS,N,'FaceColor',colors(4,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_4r,Vq,Xq,'-r','Color',colors(4,:));
axis(sfax2_4r,'off');
ylim(sfax2_4r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_4u,BINS,N,'FaceColor',colors(4,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_4u,Xq,Vq,'-r','Color',colors(4,:));
axis(sfax2_4u,'off');
xlim(sfax2_4u,[-1,1]);
title(sfax2_4u,'Working memory','FontSize',10,'FontWeight','n');

% Language
% --------
sfax2_5 = axes(sfig);
sfax2_5.Units = 'normalized';
sfax2_5.Position = [0.34,0.48,0.10,0.15];
sfax2_5r = axes(sfig);
sfax2_5r.Units = 'normalized';
sfax2_5r.Position = [0.445,0.48,0.02,0.15];
sfax2_5u = axes(sfig);
sfax2_5u.Units = 'normalized';
sfax2_5u.Position = [0.34,0.64,0.10,0.04];
hold(sfax2_5,'on');
hold(sfax2_5r,'on');
hold(sfax2_5u,'on');

x = corr_geom_beta_spin.language_math_story;
y = corr_orig_rotate_spin.language_math_story;
line(sfax2_5,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_5,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_5,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(5,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_5,[-1,1]); ylim(sfax2_5,[-1,1]);
ylabel(sfax2_5,sprintf('Correlation with\nrotated maps'))
xlabel(sfax2_5,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_5,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(5,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_5,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(5,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_5r,BINS,N,'FaceColor',colors(5,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_5r,Vq,Xq,'-r','Color',colors(5,:));
axis(sfax2_5r,'off');
ylim(sfax2_5r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_5u,BINS,N,'FaceColor',colors(5,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_5u,Xq,Vq,'-r','Color',colors(5,:));
axis(sfax2_5u,'off');
xlim(sfax2_5u,[-1,1]);
title(sfax2_5u,'Language','FontSize',10,'FontWeight','n');

% Emotion
% -------
sfax2_6 = axes(sfig);
sfax2_6.Units = 'normalized';
sfax2_6.Position = [0.50,0.48,0.10,0.15];
sfax2_6r = axes(sfig);
sfax2_6r.Units = 'normalized';
sfax2_6r.Position = [0.605,0.48,0.02,0.15];
sfax2_6u = axes(sfig);
sfax2_6u.Units = 'normalized';
sfax2_6u.Position = [0.50,0.64,0.10,0.04];
hold(sfax2_6,'on');
hold(sfax2_6r,'on');
hold(sfax2_6u,'on');

x = corr_geom_beta_spin.emotion_faces_shapes;
y = corr_orig_rotate_spin.emotion_faces_shapes;
line(sfax2_6,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_6,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_6,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(6,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_6,[-1,1]); ylim(sfax2_6,[-1,1]);
xlabel(sfax2_6,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_6,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(6,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_6,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(6,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_6r,BINS,N,'FaceColor',colors(6,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_6r,Vq,Xq,'-r','Color',colors(6,:));
axis(sfax2_6r,'off');
ylim(sfax2_6r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_6u,BINS,N,'FaceColor',colors(6,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_6u,Xq,Vq,'-r','Color',colors(6,:));
axis(sfax2_6u,'off');
xlim(sfax2_6u,[-1,1]);
title(sfax2_6u,'Emotion','FontSize',10,'FontWeight','n');

% Relational
% ----------
sfax2_7 = axes(sfig);
sfax2_7.Units = 'normalized';
sfax2_7.Position = [0.66,0.48,0.10,0.15];
sfax2_7r = axes(sfig);
sfax2_7r.Units = 'normalized';
sfax2_7r.Position = [0.765,0.48,0.02,0.15];
sfax2_7u = axes(sfig);
sfax2_7u.Units = 'normalized';
sfax2_7u.Position = [0.66,0.64,0.10,0.04];
hold(sfax2_7,'on');
hold(sfax2_7r,'on');
hold(sfax2_7u,'on');

x = corr_geom_beta_spin.relational_match_rel;
y = corr_orig_rotate_spin.relational_match_rel;
line(sfax2_7,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_7,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_7,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(7,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_7,[-1,1]); ylim(sfax2_7,[-1,1]);
xlabel(sfax2_7,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_7,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(7,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_7,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(7,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_7r,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_7r,Vq,Xq,'-r','Color',colors(7,:));
axis(sfax2_7r,'off');
ylim(sfax2_7r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_7u,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_7u,Xq,Vq,'-r','Color',colors(7,:));
axis(sfax2_7u,'off');
xlim(sfax2_7u,[-1,1]);
title(sfax2_7u,'Relational','FontSize',10,'FontWeight','n');

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure3.pdf');
%% Supplementary figure 4

% Blank figure
% ------------
sfig = figure(5);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 4');

% ---------------------------------------------------------------------------
% d: Reconstruction accuracy with rotated task activity maps (5000 rotations)
% ---------------------------------------------------------------------------
% Task labels
% -----------
str_tasks = {'emotion_faces_shapes','gambling_punish_reward','language_math_story','motor_cue_avg','relational_match_rel','social_tom_random','wm_2bk_0bk'};
str_tasks_fig = {'Emotion','Gambling','Language','Motor','Relational','Social',sprintf('Working\nmemory')};
ax4 = axes(sfig);
ax4.Units = 'normalized';
ax4.Position = [0.40,0.13,0.60,0.55];
hold on
for nTask = 1:7
    X0 = 0;
    for nMode = [1,2,5,6]
        X0 = X0 + 1;
        X1 = nTask-0.24+(X0-1)*0.16;
        val_no_rot = recon_acc_without_rotate_n200(nTask,nMode);
        eval(sprintf('val_rot    = recon_acc_n200_rot5000.%s(:,nMode);',str_tasks{nTask}));
        [X] = plot_violin_scatter(val_rot,X1,40,0.1,0);
        scatter(ax4,X,val_rot,5,'Marker','o','MarkerFaceColor',colors2(nMode,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
        hbp1 = boxplot(ax4,val_rot,'Positions',X1,'Orientation','vertical','Width',0.07,'Colors',colors2(nMode,:)*0.5,'outliersize',1.5,'symbol','k+');
        set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors2(nMode,:)*0.5);
        line([X1-0.1,X1+0.1],[1,1]*val_no_rot,'Color',colors2(3,:)*0.8,'LineWidth',1,'LineStyle','-');
        [p,h,s] = signrank(val_rot - val_no_rot);
        p = 1 - (1 - p)^7; % Bonferroni correction
        e = s.zval/sqrt(length(val_rot));
        m = max(val_rot);
        if p < 0.05 && p >=0.0001
            if median(val_rot) < val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','none');
            elseif median(val_rot) > val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','none');
            end
        elseif p < 0.0001
            if median(val_rot) < val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k');
            elseif median(val_rot) > val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k');
            end
        end
        disp(e)
    end
    
    % Task name
    text(nTask,0.79,str_tasks_fig{nTask},'FontSize',10,'FontWeight','n','HorizontalAlignment','center');
end
xlim([0.5,7.5]);
ylim([0.80,1]);
set(ax4,'FontSize',10,'FontWeight','n','XColor','w','box','off');
ax4.Box = 'off';
ylabel('Reconstruction accuracy')
scatter(ax4,3.15,0.885,15,'Marker','o','MarkerFaceColor',colors2(1,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.875,15,'Marker','o','MarkerFaceColor',colors2(2,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.865,15,'Marker','o','MarkerFaceColor',colors2(5,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.855,15,'Marker','o','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
line(ax4,[3.12,3.18],[1,1]*0.845,'Color',colors2(3,:)*0.8,'LineWidth',1,'LineStyle','-');
scatter(ax4,3.15,0.835,15,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(ax4,3.15,0.825,15,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(ax4,3.15,0.815,12,'Marker','+','MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',1.5);
text(ax4,3.2,0.885,'Geometric','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.875,'Connectome','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.865,'Parcel','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.855,'Parcel + Connectome','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.845,'Original','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.835,'Corrected \itp\rm < 0.0001 (overfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.825,'Corrected \itp\rm < 0.0001 (underfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.815,'Outliers','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
ax4.Position = [0.40,0.13,0.60,0.55];

ax4_brain1 = axes(sfig);
ax4_brain1.Units = 'normalized';
ax4_brain1.Position = [0.427,0.02,0.03,0.06];
cdata = zstat_avg_spin.emotion_faces_shapes;
patch(ax4_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain1,bluewhitered);
material(ax4_brain1,'dull');
axis(ax4_brain1,'vis3d');
axis(ax4_brain1,'off');
axis(ax4_brain1,'image');
view(ax4_brain1,270,0);
camlight(ax4_brain1,'headlight')

ax4_brain2 = axes(sfig);
ax4_brain2.Units = 'normalized';
ax4_brain2.Position = [0.513,0.02,0.03,0.06];
cdata = zstat_avg_spin.gambling_punish_reward;
patch(ax4_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain2,bluewhitered);
material(ax4_brain2,'dull');
axis(ax4_brain2,'vis3d');
axis(ax4_brain2,'off');
axis(ax4_brain2,'image');
view(ax4_brain2,270,0);
camlight(ax4_brain2,'headlight')

ax4_brain3 = axes(sfig);
ax4_brain3.Units = 'normalized';
ax4_brain3.Position = [0.599,0.02,0.03,0.06];
cdata = zstat_avg_spin.language_math_story;
patch(ax4_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain3,bluewhitered);
material(ax4_brain3,'dull');
axis(ax4_brain3,'vis3d');
axis(ax4_brain3,'off');
axis(ax4_brain3,'image');
view(ax4_brain3,270,0);
camlight(ax4_brain3,'headlight')

ax4_brain4 = axes(sfig);
ax4_brain4.Units = 'normalized';
ax4_brain4.Position = [0.685,0.02,0.03,0.06];
cdata = zstat_avg_spin.motor_cue_avg;
patch(ax4_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain4,bluewhitered);
material(ax4_brain4,'dull');
axis(ax4_brain4,'vis3d');
axis(ax4_brain4,'off');
axis(ax4_brain4,'image');
view(ax4_brain4,270,0);
camlight(ax4_brain4,'headlight')

ax4_brain5 = axes(sfig);
ax4_brain5.Units = 'normalized';
ax4_brain5.Position = [0.771,0.02,0.03,0.06];
cdata = zstat_avg_spin.relational_match_rel;
patch(ax4_brain5, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain5,bluewhitered);
material(ax4_brain5,'dull');
axis(ax4_brain5,'vis3d');
axis(ax4_brain5,'off');
axis(ax4_brain5,'image');
view(ax4_brain5,270,0);
camlight(ax4_brain5,'headlight')

ax4_brain6 = axes(sfig);
ax4_brain6.Units = 'normalized';
ax4_brain6.Position = [0.857,0.02,0.03,0.06];
cdata = zstat_avg_spin.social_tom_random;
patch(ax4_brain6, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain6,bluewhitered);
material(ax4_brain6,'dull');
axis(ax4_brain6,'vis3d');
axis(ax4_brain6,'off');
axis(ax4_brain6,'image');
view(ax4_brain6,270,0);
camlight(ax4_brain6,'headlight')

ax4_brain7 = axes(sfig);
ax4_brain7.Units = 'normalized';
ax4_brain7.Position = [0.943,0.02,0.03,0.06];
cdata = zstat_avg_spin.wm_2bk_0bk;
patch(ax4_brain7, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain7,bluewhitered);
material(ax4_brain7,'dull');
axis(ax4_brain7,'vis3d');
axis(ax4_brain7,'off');
axis(ax4_brain7,'image');
view(ax4_brain7,270,0);
camlight(ax4_brain7,'headlight')

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure4.pdf');
%% Supplementary figure 5ab

% Parameters
% ----------
nRotation_gambling = 3228;   % for randomization
nRotation_relational = 1288; % for randomization

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel','Parcel + Connectome'};

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
sfig = figure(6);clf;set(sfig,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 5ab');

% ---------------------
% a: Effect of rotation
% ---------------------
ax3 = axes(sfig);
ax3.Units = 'normalized';
ax3.Position = [0.06,0.12,0.19,0.25];
annotation(sfig, 'textbox', [0.01, 0.71, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
hold on
line(ax3,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(ax3,[0,0],[-0.4,0.8],'Color','k','LineStyle',':','LineWidth',1);
scatter(ax3,corr_geom_beta_moran.gambling_punish_reward,corr_orig_rotate_moran.gambling_punish_reward,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',0.3)
scatter(ax3,corr_geom_beta_moran.relational_match_rel,corr_orig_rotate_moran.relational_match_rel,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(7,:),'MarkerFaceAlpha',0.3)
xlim(ax3,[-1,1]);
ylim(ax3,[-0.4,0.8]);
xlabel(ax3,'Correlation between \beta values')
ylabel(ax3,'Correlation with rotated maps')

x = corr_geom_beta_moran.gambling_punish_reward;
y = corr_orig_rotate_moran.gambling_punish_reward;
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(ax3,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(3,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(ax3,-0.95,0.74,sprintf('Pearson r = %.03f',r),'Color',colors(3,:),'HorizontalAlignment','left');
r = corr(x',y','type','Spearman');
text(ax3,-0.95,0.64,sprintf('Spearman r = %.03f',r),'Color',colors(3,:),'HorizontalAlignment','left');

x = corr_geom_beta_moran.relational_match_rel;
y = corr_orig_rotate_moran.relational_match_rel;
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(ax3,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(7,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(ax3,-0.95,0.51,sprintf('Pearson r = %.03f',r),'Color',colors(7,:),'HorizontalAlignment','left');
r = corr(x',y','type','Spearman');
text(ax3,-0.95,0.41,sprintf('Spearman r = %.03f',r),'Color',colors(7,:),'HorizontalAlignment','left');

ax3_r = axes(sfig);
ax3_r.Units = 'normalized';
ax3_r.Position = [0.255,0.12,0.04,0.25];
Y = corr_orig_rotate_moran.relational_match_rel;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(ax3_r,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_r,Vq,Xq,'-r','Color',colors(7,:));
axis(ax3_r,'off');
ylim(ax3_r,[-0.4,0.8]);

Y = corr_orig_rotate_moran.gambling_punish_reward;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(ax3_r,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_r,Vq,Xq,'-k','Color',colors(3,:))
axis(ax3_r,'off');
ylim(ax3_r,[-0.4,0.8]);

ax3_u = axes(sfig);
ax3_u.Units = 'normalized';
ax3_u.Position = [0.06,0.38,0.19,0.08];
X = corr_geom_beta_moran.relational_match_rel;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(ax3_u,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_u,Xq,Vq,'-r','Color',colors(7,:));
axis(ax3_u,'off');
xlim(ax3_u,[-1,1]);

X = corr_geom_beta_moran.gambling_punish_reward;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(ax3_u,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
hold on;
plot(ax3_u,Xq,Vq,'-k','Color',colors(3,:));
axis(ax3_u,'off');
xlim(ax3_u,[-1,1]);

ax3_text1 = axes(sfig);
ax3_text1.Units = 'normalized';
ax3_text1.Position = [0.02,0.62,0.01,0.01];
xlim(ax3_text1,[-1,1]);ylim(ax3_text1,[-1,1]);axis(ax3_text1,'off');
text(ax3_text1,0,0,'Gambling','FontSize',10,'FontWeight','n','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','Color',colors(3,:));

ax3_text11 = axes(sfig);
ax3_text11.Units = 'normalized';
ax3_text11.Position = [0.125,0.67,0.01,0.01];
xlim(ax3_text11,[-1,1]);ylim(ax3_text11,[-1,1]);axis(ax3_text11,'off');
text(ax3_text11,0,0,'Original','FontSize',10,'FontWeight','n','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','Color','k');

ax3_text12 = axes(sfig);
ax3_text12.Units = 'normalized';
ax3_text12.Position = [0.245,0.67,0.01,0.01];
xlim(ax3_text12,[-1,1]);ylim(ax3_text12,[-1,1]);axis(ax3_text12,'off');
text(ax3_text12,0,0,'Randomized map','FontSize',10,'FontWeight','n','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','Color','k');

ax3_text2 = axes(sfig);
ax3_text2.Units = 'normalized';
ax3_text2.Position = [0.02,0.52,0.01,0.01];
xlim(ax3_text2,[-1,1]);ylim(ax3_text2,[-1,1]);axis(ax3_text2,'off');
text(ax3_text2,0,0,'Relational','FontSize',10,'FontWeight','n','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','Color',colors(7,:));

ax3_brain1 = axes(sfig);
ax3_brain1.Units = 'normalized';
cdata = zstat_avg_moran.gambling_punish_reward; cdata(isnan(cdata)) = 0;
patch(ax3_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain1,bluewhitered);
material(ax3_brain1,'dull');
axis(ax3_brain1,'vis3d');
axis(ax3_brain1,'off');
axis(ax3_brain1,'image');
view(ax3_brain1,270,0);
camlight(ax3_brain1,'headlight')
ax3_brain1.Position = [0.08,0.58,0.04,0.08];

ax3_brain1_flip = axes(sfig);
ax3_brain1_flip.Units = 'normalized';
cdata = zstat_avg_moran.gambling_punish_reward; cdata(isnan(cdata)) = 0;
patch(ax3_brain1_flip, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain1_flip,bluewhitered);
material(ax3_brain1_flip,'dull');
axis(ax3_brain1_flip,'vis3d');
axis(ax3_brain1_flip,'off');
axis(ax3_brain1_flip,'image');
view(ax3_brain1_flip,90,0);
camlight(ax3_brain1_flip,'headlight')
ax3_brain1_flip.Position = [0.14,0.58,0.04,0.08];

ax3_brain2 = axes(sfig);
ax3_brain2.Units = 'normalized';
cdata = squeeze(zstat_random.gambling_punish_reward{1,1}(:,1,nRotation_gambling));
patch(ax3_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain2,bluewhitered);
material(ax3_brain2,'dull')
axis(ax3_brain2,'vis3d')
axis(ax3_brain2,'off');
axis(ax3_brain2,'image');
view(ax3_brain2,270,0);
camlight(ax3_brain2,'headlight')
ax3_brain2.Position = [0.20,0.58,0.04,0.08];

ax3_brain2_flip = axes(sfig);
ax3_brain2_flip.Units = 'normalized';
cdata = squeeze(zstat_random.gambling_punish_reward{1,1}(:,1,nRotation_gambling));
patch(ax3_brain2_flip, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain2_flip,bluewhitered);
material(ax3_brain2_flip,'dull')
axis(ax3_brain2_flip,'vis3d')
axis(ax3_brain2_flip,'off');
axis(ax3_brain2_flip,'image');
view(ax3_brain2_flip,90,0);
camlight(ax3_brain2_flip,'headlight')
ax3_brain2_flip.Position = [0.26,0.58,0.04,0.08];

scatter(ax3,corr_geom_beta_moran.gambling_punish_reward(1,nRotation_gambling),corr_orig_rotate_moran.gambling_punish_reward(1,nRotation_gambling),20,'Marker','o','MarkerEdgeColor',colors(3,:)*0.6,'LineWidth',2);
ax3.Position = [0.06,0.12,0.20,0.26];

ax3_brain3 = axes(sfig);
ax3_brain3.Units = 'normalized';
cdata = zstat_avg_moran.relational_match_rel; cdata(isnan(cdata)) = 0;
patch(ax3_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain3,bluewhitered);
material(ax3_brain3,'dull');
axis(ax3_brain3,'vis3d');
axis(ax3_brain3,'off');
axis(ax3_brain3,'image');
view(ax3_brain3,270,0);
camlight(ax3_brain3,'headlight')
ax3_brain3.Position = [0.08,0.48,0.04,0.08];

ax3_brain3_flip = axes(sfig);
ax3_brain3_flip.Units = 'normalized';
cdata = zstat_avg_moran.relational_match_rel; cdata(isnan(cdata)) = 0;
patch(ax3_brain3_flip, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain3_flip,bluewhitered);
material(ax3_brain3_flip,'dull');
axis(ax3_brain3_flip,'vis3d');
axis(ax3_brain3_flip,'off');
axis(ax3_brain3_flip,'image');
view(ax3_brain3_flip,90,0);
camlight(ax3_brain3_flip,'headlight')
ax3_brain3_flip.Position = [0.14,0.48,0.04,0.08];

ax3_brain4 = axes(sfig);
ax3_brain4.Units = 'normalized';
cdata = squeeze(zstat_random.relational_match_rel{1,1}(:,1,nRotation_relational));
patch(ax3_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain4,bluewhitered);
material(ax3_brain4,'dull')
axis(ax3_brain4,'vis3d')
axis(ax3_brain4,'off');
axis(ax3_brain4,'image');
view(ax3_brain4,270,0);
camlight(ax3_brain4,'headlight')
ax3_brain4.Position = [0.20,0.48,0.04,0.08];

ax3_brain4_flip = axes(sfig);
ax3_brain4_flip.Units = 'normalized';
cdata = squeeze(zstat_random.relational_match_rel{1,1}(:,1,nRotation_relational));
patch(ax3_brain4_flip, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax3_brain4_flip,bluewhitered);
material(ax3_brain4_flip,'dull')
axis(ax3_brain4_flip,'vis3d')
axis(ax3_brain4_flip,'off');
axis(ax3_brain4_flip,'image');
view(ax3_brain4_flip,90,0);
camlight(ax3_brain4_flip,'headlight')
ax3_brain4_flip.Position = [0.26,0.48,0.04,0.08];

scatter(ax3,corr_geom_beta_moran.relational_match_rel(1,nRotation_relational),corr_orig_rotate_moran.relational_match_rel(1,nRotation_relational),20,'Marker','o','MarkerEdgeColor',colors(7,:)*0.6,'LineWidth',2);
ax3.Position = [0.06,0.12,0.20,0.26];

ax3_accu = axes(sfig);
ax3_accu.Units = 'normalized';
hold on
Y1 = recon_acc_n200_rnd5000_no_add.gambling_punish_reward(:,1);
Y2 = recon_acc_n200_rnd5000_no_add.relational_match_rel(:,1);
val_rot_no_add = Y1;
[X] = plot_violin_scatter(val_rot_no_add,1,40,0.5,0);
scatter(ax3_accu,X,val_rot_no_add,5,'Marker','o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
hbp1 = boxplot(ax3_accu,val_rot_no_add,'Positions',1,'Orientation','vertical','Width',0.3,'Colors',colors(3,:)*0.5,'outliersize',1.5,'symbol','k+');
set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors(3,:)*0.5);
val_rot_no_add = Y2;
[X] = plot_violin_scatter(val_rot_no_add,2,40,0.5,0);
scatter(ax3_accu,X,val_rot_no_add,5,'Marker','o','MarkerFaceColor',colors(7,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
hbp1 = boxplot(ax3_accu,val_rot_no_add,'Positions',2,'Orientation','vertical','Width',0.3,'Colors',colors(7,:)*0.5,'outliersize',1.5,'symbol','k+');
set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors(7,:)*0.5);
xlim([0.5,2.5]);
ylim([0.9,1]);
axis(ax3_accu,'on');
set(ax3_accu,'FontSize',10,'FontWeight','n','XColor','w','box','off','YTick',[0.9,1.0]);
ylabel(ax3_accu,sprintf('Recon. acc.'));
text(ax3_accu,1.5,0.88,'Geometric','HorizontalAlignment','center','FontSize',10,'FontWeight','n');
line(ax3_accu,[1,1,2,2],[0.994,1,1,0.997],'Color','k');
ax3_accu.Position = [0.29,0.31,0.03,0.15];

ax3.Position = [0.06,0.12,0.19,0.25];

% ---------------------------------------------------------------------------
% b: Reconstruction accuracy with rotated task activity maps (5000 rotations)
% ---------------------------------------------------------------------------
% Task labels
% -----------
annotation(sfig, 'textbox', [0.36, 0.71, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
str_tasks = {'emotion_faces_shapes','gambling_punish_reward','language_math_story','motor_cue_avg','relational_match_rel','social_tom_random','wm_2bk_0bk'};
str_tasks_fig = {'Emotion','Gambling','Language','Motor','Relational','Social',sprintf('Working\nmemory')};
ax4 = axes(sfig);
ax4.Units = 'normalized';
ax4.Position = [0.40,0.13,0.60,0.55];
hold on
for nTask = 1:7
    X0 = 0;
    for nMode = [1,2,5,6]
        X0 = X0 + 1;
        X1 = nTask-0.24+(X0-1)*0.16;
        val_no_rot = recon_acc_without_rotate_n200(nTask,nMode);
        eval(sprintf('val_rot    = recon_acc_n200_rnd5000.%s(:,nMode);',str_tasks{nTask}));
        [X] = plot_violin_scatter(val_rot,X1,40,0.1,0);
        scatter(ax4,X,val_rot,5,'Marker','o','MarkerFaceColor',colors2(nMode,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
        hbp1 = boxplot(ax4,val_rot,'Positions',X1,'Orientation','vertical','Width',0.07,'Colors',colors2(nMode,:)*0.5,'outliersize',1.5,'symbol','k+');
        set(hbp1,{'linew'},{1.5});set(hbp1(6),'Color',colors2(nMode,:)*0.5);
        line([X1-0.1,X1+0.1],[1,1]*val_no_rot,'Color',colors2(3,:)*0.8,'LineWidth',1,'LineStyle','-');
        [p,h,s] = signrank(val_rot - val_no_rot);
        p = 1 - (1 - p)^7; % Bonferroni correction
        e = s.zval/sqrt(length(val_rot));
        m = max(val_rot);
        if p < 0.05 && p >=0.0001
            if median(val_rot) < val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','none');
            elseif median(val_rot) > val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','none');
            end
        elseif p < 0.0001
            if median(val_rot) < val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k');
            elseif median(val_rot) > val_no_rot
                scatter(ax4,X1,m + 0.004,20,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k');
            end
        end
        disp(e)
    end
    
    % Task name
    text(nTask,0.79,str_tasks_fig{nTask},'FontSize',10,'FontWeight','n','HorizontalAlignment','center');
end
xlim([0.5,7.5]);
ylim([0.80,1]);
set(ax4,'FontSize',10,'FontWeight','n','XColor','w','box','off');
ax4.Box = 'off';
ylabel('Reconstruction accuracy')
scatter(ax4,3.15,0.885,15,'Marker','o','MarkerFaceColor',colors2(1,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.875,15,'Marker','o','MarkerFaceColor',colors2(2,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.865,15,'Marker','o','MarkerFaceColor',colors2(5,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
scatter(ax4,3.15,0.855,15,'Marker','o','MarkerFaceColor',colors2(6,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
line(ax4,[3.12,3.18],[1,1]*0.845,'Color',colors2(3,:)*0.8,'LineWidth',1,'LineStyle','-');
scatter(ax4,3.15,0.835,15,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(ax4,3.15,0.825,15,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(ax4,3.15,0.815,12,'Marker','^','MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(ax4,3.15,0.805,12,'Marker','+','MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',1.5);
text(ax4,3.2,0.885,'Geometric','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.875,'Connectome','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.865,'Parcel','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.855,'Parcel + Connectome','FontSize',10,'FontWeight','n','HorizontalAlignment','left');
text(ax4,3.2,0.845,'Original','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.835,'Corrected \itp\rm < 0.0001 (overfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.825,'Corrected \itp\rm < 0.0001 (underfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.815,'Corrected \itp\rm < 0.05 (overfit)','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
text(ax4,3.2,0.805,'Outliers','HorizontalAlignment','left','FontSize',10,'FontWeight','n')
ax4.Position = [0.40,0.13,0.60,0.55];

ax4_brain1 = axes(sfig);
ax4_brain1.Units = 'normalized';
ax4_brain1.Position = [0.427,0.02,0.03,0.06];
cdata = zstat_avg_moran.emotion_faces_shapes;
patch(ax4_brain1, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain1,bluewhitered);
material(ax4_brain1,'dull');
axis(ax4_brain1,'vis3d');
axis(ax4_brain1,'off');
axis(ax4_brain1,'image');
view(ax4_brain1,270,0);
camlight(ax4_brain1,'headlight')

ax4_brain2 = axes(sfig);
ax4_brain2.Units = 'normalized';
ax4_brain2.Position = [0.513,0.02,0.03,0.06];
cdata = zstat_avg_moran.gambling_punish_reward;
patch(ax4_brain2, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain2,bluewhitered);
material(ax4_brain2,'dull');
axis(ax4_brain2,'vis3d');
axis(ax4_brain2,'off');
axis(ax4_brain2,'image');
view(ax4_brain2,270,0);
camlight(ax4_brain2,'headlight')

ax4_brain3 = axes(sfig);
ax4_brain3.Units = 'normalized';
ax4_brain3.Position = [0.599,0.02,0.03,0.06];
cdata = zstat_avg_moran.language_math_story;
patch(ax4_brain3, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain3,bluewhitered);
material(ax4_brain3,'dull');
axis(ax4_brain3,'vis3d');
axis(ax4_brain3,'off');
axis(ax4_brain3,'image');
view(ax4_brain3,270,0);
camlight(ax4_brain3,'headlight')

ax4_brain4 = axes(sfig);
ax4_brain4.Units = 'normalized';
ax4_brain4.Position = [0.685,0.02,0.03,0.06];
cdata = zstat_avg_moran.motor_cue_avg;
patch(ax4_brain4, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain4,bluewhitered);
material(ax4_brain4,'dull');
axis(ax4_brain4,'vis3d');
axis(ax4_brain4,'off');
axis(ax4_brain4,'image');
view(ax4_brain4,270,0);
camlight(ax4_brain4,'headlight')

ax4_brain5 = axes(sfig);
ax4_brain5.Units = 'normalized';
ax4_brain5.Position = [0.771,0.02,0.03,0.06];
cdata = zstat_avg_moran.relational_match_rel;
patch(ax4_brain5, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain5,bluewhitered);
material(ax4_brain5,'dull');
axis(ax4_brain5,'vis3d');
axis(ax4_brain5,'off');
axis(ax4_brain5,'image');
view(ax4_brain5,270,0);
camlight(ax4_brain5,'headlight')

ax4_brain6 = axes(sfig);
ax4_brain6.Units = 'normalized';
ax4_brain6.Position = [0.857,0.02,0.03,0.06];
cdata = zstat_avg_moran.social_tom_random;
patch(ax4_brain6, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain6,bluewhitered);
material(ax4_brain6,'dull');
axis(ax4_brain6,'vis3d');
axis(ax4_brain6,'off');
axis(ax4_brain6,'image');
view(ax4_brain6,270,0);
camlight(ax4_brain6,'headlight')

ax4_brain7 = axes(sfig);
ax4_brain7.Units = 'normalized';
ax4_brain7.Position = [0.943,0.02,0.03,0.06];
cdata = zstat_avg_moran.wm_2bk_0bk;
patch(ax4_brain7, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(ax4_brain7,bluewhitered);
material(ax4_brain7,'dull');
axis(ax4_brain7,'vis3d');
axis(ax4_brain7,'off');
axis(ax4_brain7,'image');
view(ax4_brain7,270,0);
camlight(ax4_brain7,'headlight')

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure5ab.pdf');
%% Supplementary figure 5c

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel','Parcel + Connectome'};

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
sfig = figure(7);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 5c');
annotation(sfig, 'textbox', [0.29, 0.99, 0.01, 0.01], 'string', 'c', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

% Social
% ------
sfax2_1 = axes(sfig);
sfax2_1.Units = 'normalized';
sfax2_1.Position = [0.34,0.75,0.10,0.15];
sfax2_1r = axes(sfig);
sfax2_1r.Units = 'normalized';
sfax2_1r.Position = [0.445,0.75,0.02,0.15];
sfax2_1u = axes(sfig);
sfax2_1u.Units = 'normalized';
sfax2_1u.Position = [0.34,0.91,0.10,0.04];
hold(sfax2_1,'on');
hold(sfax2_1r,'on');
hold(sfax2_1u,'on');

line(sfax2_1,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_1,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_1,corr_geom_beta_moran.social_tom_random,corr_orig_rotate_moran.social_tom_random,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(1,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_1,[-1,1]); ylim(sfax2_1,[-1,1]);
ylabel(sfax2_1,sprintf('Correlation with\nrotated maps'))
x = corr_geom_beta_moran.social_tom_random;
y = corr_orig_rotate_moran.social_tom_random;
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_1,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(1,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_1,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(1,:),'HorizontalAlignment','left');

Y = corr_orig_rotate_moran.social_tom_random;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_1r,BINS,N,'FaceColor',colors(1,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_1r,Vq,Xq,'-r','Color',colors(1,:));
axis(sfax2_1r,'off');
ylim(sfax2_1r,[-1,1]);

X = corr_geom_beta_moran.social_tom_random;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_1u,BINS,N,'FaceColor',colors(1,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_1u,Xq,Vq,'-r','Color',colors(1,:));
axis(sfax2_1u,'off');
xlim(sfax2_1u,[-1,1]);
title(sfax2_1u,'Social','FontSize',10,'FontWeight','n');

% Motor
% -----
sfax2_2 = axes(sfig);
sfax2_2.Units = 'normalized';
sfax2_2.Position = [0.50,0.75,0.10,0.15];
sfax2_2r = axes(sfig);
sfax2_2r.Units = 'normalized';
sfax2_2r.Position = [0.605,0.75,0.02,0.15];
sfax2_2u = axes(sfig);
sfax2_2u.Units = 'normalized';
sfax2_2u.Position = [0.50,0.91,0.10,0.04];
hold(sfax2_2,'on');
hold(sfax2_2r,'on');
hold(sfax2_2u,'on');

line(sfax2_2,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_2,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
x = corr_geom_beta_moran.motor_cue_avg;
y = corr_orig_rotate_moran.motor_cue_avg;
scatter(sfax2_2,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_2,[-1,1]); ylim(sfax2_2,[-1,1]);
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_2,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(2,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_2,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(2,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_2r,BINS,N,'FaceColor',colors(2,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_2r,Vq,Xq,'-r','Color',colors(2,:));
axis(sfax2_2r,'off');
ylim(sfax2_2r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_2u,BINS,N,'FaceColor',colors(2,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_2u,Xq,Vq,'-r','Color',colors(2,:));
axis(sfax2_2u,'off');
xlim(sfax2_2u,[-1,1]);
title(sfax2_2u,'Motor','FontSize',10,'FontWeight','n');

% Gambling
% --------
sfax2_3 = axes(sfig);
sfax2_3.Units = 'normalized';
sfax2_3.Position = [0.66,0.75,0.10,0.15];
sfax2_3r = axes(sfig);
sfax2_3r.Units = 'normalized';
sfax2_3r.Position = [0.765,0.75,0.02,0.15];
sfax2_3u = axes(sfig);
sfax2_3u.Units = 'normalized';
sfax2_3u.Position = [0.66,0.91,0.10,0.04];
hold(sfax2_3,'on');
hold(sfax2_3r,'on');
hold(sfax2_3u,'on');

line(sfax2_3,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_3,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
x = corr_geom_beta_moran.gambling_punish_reward;
y = corr_orig_rotate_moran.gambling_punish_reward;
scatter(sfax2_3,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_3,[-1,1]); ylim(sfax2_3,[-1,1]);
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_3,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(3,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_3,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(3,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_3r,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_3r,Vq,Xq,'-r','Color',colors(3,:));
axis(sfax2_3r,'off');
ylim(sfax2_3r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_3u,BINS,N,'FaceColor',colors(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_3u,Xq,Vq,'-r','Color',colors(3,:));
axis(sfax2_3u,'off');
xlim(sfax2_3u,[-1,1]);
title(sfax2_3u,'Gambling','FontSize',10,'FontWeight','n');

% Working memory
% --------------
sfax2_4 = axes(sfig);
sfax2_4.Units = 'normalized';
sfax2_4.Position = [0.82,0.75,0.10,0.15];
sfax2_4r = axes(sfig);
sfax2_4r.Units = 'normalized';
sfax2_4r.Position = [0.925,0.75,0.02,0.15];
sfax2_4u = axes(sfig);
sfax2_4u.Units = 'normalized';
sfax2_4u.Position = [0.82,0.91,0.10,0.04];
hold(sfax2_4,'on');
hold(sfax2_4r,'on');
hold(sfax2_4u,'on');

line(sfax2_4,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_4,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
x = corr_geom_beta_moran.wm_2bk_0bk;
y = corr_orig_rotate_moran.wm_2bk_0bk;
scatter(sfax2_4,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(4,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_4,[-1,1]); ylim(sfax2_4,[-1,1]);
xlabel(sfax2_4,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_4,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(4,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_4,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(4,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_4r,BINS,N,'FaceColor',colors(4,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_4r,Vq,Xq,'-r','Color',colors(4,:));
axis(sfax2_4r,'off');
ylim(sfax2_4r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_4u,BINS,N,'FaceColor',colors(4,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_4u,Xq,Vq,'-r','Color',colors(4,:));
axis(sfax2_4u,'off');
xlim(sfax2_4u,[-1,1]);
title(sfax2_4u,'Working memory','FontSize',10,'FontWeight','n');

% Language
% --------
sfax2_5 = axes(sfig);
sfax2_5.Units = 'normalized';
sfax2_5.Position = [0.34,0.48,0.10,0.15];
sfax2_5r = axes(sfig);
sfax2_5r.Units = 'normalized';
sfax2_5r.Position = [0.445,0.48,0.02,0.15];
sfax2_5u = axes(sfig);
sfax2_5u.Units = 'normalized';
sfax2_5u.Position = [0.34,0.64,0.10,0.04];
hold(sfax2_5,'on');
hold(sfax2_5r,'on');
hold(sfax2_5u,'on');

x = corr_geom_beta_moran.language_math_story;
y = corr_orig_rotate_moran.language_math_story;
line(sfax2_5,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_5,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_5,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(5,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_5,[-1,1]); ylim(sfax2_5,[-1,1]);
ylabel(sfax2_5,sprintf('Correlation with\nrotated maps'))
xlabel(sfax2_5,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_5,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(5,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_5,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(5,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_5r,BINS,N,'FaceColor',colors(5,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_5r,Vq,Xq,'-r','Color',colors(5,:));
axis(sfax2_5r,'off');
ylim(sfax2_5r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_5u,BINS,N,'FaceColor',colors(5,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_5u,Xq,Vq,'-r','Color',colors(5,:));
axis(sfax2_5u,'off');
xlim(sfax2_5u,[-1,1]);
title(sfax2_5u,'Language','FontSize',10,'FontWeight','n');

% Emotion
% -------
sfax2_6 = axes(sfig);
sfax2_6.Units = 'normalized';
sfax2_6.Position = [0.50,0.48,0.10,0.15];
sfax2_6r = axes(sfig);
sfax2_6r.Units = 'normalized';
sfax2_6r.Position = [0.605,0.48,0.02,0.15];
sfax2_6u = axes(sfig);
sfax2_6u.Units = 'normalized';
sfax2_6u.Position = [0.50,0.64,0.10,0.04];
hold(sfax2_6,'on');
hold(sfax2_6r,'on');
hold(sfax2_6u,'on');

x = corr_geom_beta_moran.emotion_faces_shapes;
y = corr_orig_rotate_moran.emotion_faces_shapes;
line(sfax2_6,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_6,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_6,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(6,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_6,[-1,1]); ylim(sfax2_6,[-1,1]);
xlabel(sfax2_6,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_6,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(6,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_6,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(6,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_6r,BINS,N,'FaceColor',colors(6,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_6r,Vq,Xq,'-r','Color',colors(6,:));
axis(sfax2_6r,'off');
ylim(sfax2_6r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_6u,BINS,N,'FaceColor',colors(6,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_6u,Xq,Vq,'-r','Color',colors(6,:));
axis(sfax2_6u,'off');
xlim(sfax2_6u,[-1,1]);
title(sfax2_6u,'Emotion','FontSize',10,'FontWeight','n');

% Relational
% ----------
sfax2_7 = axes(sfig);
sfax2_7.Units = 'normalized';
sfax2_7.Position = [0.66,0.48,0.10,0.15];
sfax2_7r = axes(sfig);
sfax2_7r.Units = 'normalized';
sfax2_7r.Position = [0.765,0.48,0.02,0.15];
sfax2_7u = axes(sfig);
sfax2_7u.Units = 'normalized';
sfax2_7u.Position = [0.66,0.64,0.10,0.04];
hold(sfax2_7,'on');
hold(sfax2_7r,'on');
hold(sfax2_7u,'on');

x = corr_geom_beta_moran.relational_match_rel;
y = corr_orig_rotate_moran.relational_match_rel;
line(sfax2_7,[-1,1],[0,0],'Color','k','LineStyle',':','LineWidth',1);
line(sfax2_7,[0,0],[-1,1],'Color','k','LineStyle',':','LineWidth',1);
scatter(sfax2_7,x,y,5,'Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(7,:),'MarkerFaceAlpha',0.3)
xlim(sfax2_7,[-1,1]); ylim(sfax2_7,[-1,1]);
xlabel(sfax2_7,'Correlation between \beta values')
lm = fitlm(x, y);
B = lm.Coefficients.Estimate;
line(sfax2_7,[-1,1],B(1)+B(2)*([-1,1]),'Color',colors(7,:)*0.8,'LineWidth',1);
r = corr(x',y','type','Pearson');
text(sfax2_7,-0.95,0.8,sprintf('r = %.03f',r),'Color',colors(7,:),'HorizontalAlignment','left');

Y = y;
[N,EDGES] = histcounts(Y,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = barh(sfax2_7r,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_7r,Vq,Xq,'-r','Color',colors(7,:));
axis(sfax2_7r,'off');
ylim(sfax2_7r,[-1,1]);

X = x;
[N,EDGES] = histcounts(X,-1:0.1/4:1);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_7u,BINS,N,'FaceColor',colors(7,:)*0.5,'EdgeColor','none','FaceAlpha',0.4);
Xq = BINS;
Vq = movmean(N,5);
plot(sfax2_7u,Xq,Vq,'-r','Color',colors(7,:));
axis(sfax2_7u,'off');
xlim(sfax2_7u,[-1,1]);
title(sfax2_7u,'Relational','FontSize',10,'FontWeight','n');

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure5c.pdf');
%% Supplementary figure 6

% Mode labels
% -----------
str_modes = {'Geometric','Connectome',sprintf('Connectome (density matched)'),'EDR','Parcel','Parcel + Connectome'};

% Calculate binscatter
% --------------------
eigenmodes_names = fields(l2_err);
xbin = 0:0.2/50:0.2;
ybin = 0:5/50:5.0;
bin_maps = [];
for nMode = 1:length(eigenmodes_names)
    emode = eigenmodes_names{nMode};
    X = l2_err.(emode)(:,1);
    Y = abs(l2_err.(emode)(:,2));
    temp_map = zeros(numel(xbin)-1,numel(ybin)-1);
    for nx = 1:numel(xbin)-1
        lgc_x = X >= xbin(nx) & X < xbin(nx+1);
        for ny = 1:numel(ybin)-1
            lgc_y = Y(lgc_x) >= ybin(ny) & Y(lgc_x) < ybin(ny+1);
            temp_map(nx,ny) = sum(lgc_y);
        end
    end
    bin_maps.(emode) = temp_map;
end

% Colormap
% --------
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0.5,0.5,0.5; cmap1([4,5,1:3,6:7],:)];

% Blank figure
% ------------
sfig = figure(8);clf;set(gcf,'Color','w','Position',[1,1,1000,500],'Name','Suppl. Figure 6');

% ---------------------------
% a-d: L2-norm vs. Abs. error
% ---------------------------
bin_size = 0.2/400;
cmap = pink(500); cmap = [1,1,1;cmap(end-50:-1:1,:)];
sfax1 = axes(sfig);
sfax1.Units = 'normalized';
sfax1.Position = [0.04,0.65,0.15,0.3];
annotation(sfig, 'textbox', [0.01, 0.99, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
temp_map = bin_maps.Geometric';
imagesc(sfax1,temp_map);colormap(sfax1,cmap);axis(sfax1,'image','xy');
hcb = colorbar(sfax1);
sfax1.Position = [0.04,0.65,0.15,0.3];
set(sfax1,'XTick',1,'XTickLabel','')
set(sfax1,'YTick',(11:10:length(ybin))-0.5,'YTickLabel',ybin(11:10:length(ybin)))
ylabel(sfax1,'Absolute error');
set(sfax1,'FontSize',10,'FontWeight','n');
title('Geometric','FontSize',10,'FontWeight','n')
cr = corr(l2_err.Geometric(:,1), abs(l2_err.Geometric(:,2)), 'Rows', 'complete');
text(sfax1,2,48,sprintf('r = %.03f',cr),'FontSize',10,'FontWeight','n');

sfax1_brain = axes(sfig);
sfax1_brain.Units = 'normalized';
sfax1_brain.Position = [0.11,0.81,0.06,0.12];
cdata = eigenmodes_l2norm.Geometric;
patch(sfax1_brain, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(sfax1_brain,cmap);
caxis(sfax1_brain,[0.02,0.15]);
material(sfax1_brain,'dull');
axis(sfax1_brain,'vis3d');
axis(sfax1_brain,'off');
axis(sfax1_brain,'image');
view(sfax1_brain,270,0);
camlight(sfax1_brain,'headlight')

sfax1_hist = axes(sfig);
sfax1_hist.Units = 'normalized';
sfax1_hist.Position = [0.04,0.58,0.15,0.06];
Y = eigenmodes_l2norm.Geometric;
[N,EDGES] = histcounts(Y,0:bin_size:0.2);
BINS = EDGES(2:end)-(EDGES(2) - EDGES(1))/2;
h = bar(sfax1_hist,BINS,N,'FaceColor',colors2(1,:)*0.5,'EdgeColor','none','FaceAlpha',0.3);
Xq = BINS;
Vq = movmean(N,5);
hold(sfax1_hist,'on');
plot(sfax1_hist,Xq,Vq,'-r','Color',colors2(1,:),'LineWidth',1);
set(sfax1_hist,'box','off','YColor','w','box','on');
ylim(sfax1_hist,[0,4000])
xlabel(sfax1_hist,'L2-norm');

sfax2 = axes(sfig);
sfax2.Units = 'normalized';
sfax2.Position = [0.28,0.65,0.15,0.3];
annotation(sfig, 'textbox', [0.25, 0.99, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
temp_map = bin_maps.Connectome';
imagesc(sfax2,temp_map);colormap(sfax2,cmap);axis(sfax2,'image','xy');
hcb = colorbar(sfax2);
sfax2.Position = [0.28,0.65,0.15,0.3];
set(sfax2,'XTick',1,'XTickLabel','')
set(sfax2,'YTick',(11:10:length(ybin))-0.5,'YTickLabel',ybin(11:10:length(ybin)))
xlabel('L2-norm');
set(sfax2,'FontSize',10,'FontWeight','n')
title('Connectome','FontSize',10,'FontWeight','n')
cr = corr(l2_err.Connectome(:,1), abs(l2_err.Connectome(:,2)), 'Rows', 'complete');
text(sfax2,2,48,sprintf('r = %.03f',cr),'FontSize',10,'FontWeight','n');

sfax2_brain = axes(sfig);
sfax2_brain.Units = 'normalized';
sfax2_brain.Position = [0.35,0.81,0.06,0.12];
cdata = eigenmodes_l2norm.Connectome;
patch(sfax2_brain, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(sfax2_brain,cmap);
caxis(sfax2_brain,[0.02,0.15]);
material(sfax2_brain,'dull');
axis(sfax2_brain,'vis3d');
axis(sfax2_brain,'off');
axis(sfax2_brain,'image');
view(sfax2_brain,270,0);
camlight(sfax2_brain,'headlight')

sfax2_hist = axes(sfig);
sfax2_hist.Units = 'normalized';
sfax2_hist.Position = [0.28,0.58,0.15,0.06];
Y = eigenmodes_l2norm.Connectome;
[N,EDGES] = histcounts(Y,0:bin_size:0.2);
BINS = EDGES(2:end) -(EDGES(2) - EDGES(1))/2;
h = bar(sfax2_hist,BINS,N,'FaceColor',colors2(2,:)*0.5,'EdgeColor','none','FaceAlpha',0.3);
Xq = BINS;
Vq = movmean(N,5);
hold(sfax2_hist,'on');
plot(sfax2_hist,Xq,Vq,'-r','Color',colors2(2,:),'LineWidth',1);
set(sfax2_hist,'box','off','YColor','w','box','on');
ylim(sfax2_hist,[0,350])
xlabel(sfax2_hist,'L2-norm');

sfax3 = axes(sfig);
sfax3.Units = 'normalized';
sfax3.Position = [0.53,0.65,0.15,0.3];
annotation(sfig, 'textbox', [0.50, 0.99, 0.01, 0.01], 'string', 'c', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
temp_map = bin_maps.Connectome_DM';
imagesc(sfax3,temp_map);colormap(sfax3,cmap);axis(sfax3,'image','xy');
hcb = colorbar(sfax3);
sfax3.Position = [0.53,0.65,0.15,0.3];
set(sfax3,'XTick',1,'XTickLabel','')
set(sfax3,'YTick',(11:10:length(ybin))-0.5,'YTickLabel',ybin(11:10:length(ybin)))
xlabel('L2-norm');
set(sfax3,'FontSize',10,'FontWeight','n')
title('Connectome (density matched)','FontSize',10,'FontWeight','n')
cr = corr(l2_err.Connectome_DM(:,1), abs(l2_err.Connectome_DM(:,2)), 'Rows', 'complete');
text(sfax3,2,48,sprintf('r = %.03f',cr),'FontSize',10,'FontWeight','n');

sfax3_brain = axes(sfig);
sfax3_brain.Units = 'normalized';
sfax3_brain.Position = [0.60,0.81,0.06,0.12];
cdata = eigenmodes_l2norm.Connectome_DM;
patch(sfax3_brain, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(sfax3_brain,cmap);
caxis(sfax3_brain,[0.02,0.15]);
material(sfax3_brain,'dull');
axis(sfax3_brain,'vis3d');
axis(sfax3_brain,'off');
axis(sfax3_brain,'image');
view(sfax3_brain,270,0);
camlight(sfax3_brain,'headlight')

sfax3_hist = axes(sfig);
sfax3_hist.Units = 'normalized';
sfax3_hist.Position = [0.53,0.58,0.15,0.06];
Y = eigenmodes_l2norm.Connectome_DM;
[N,EDGES] = histcounts(Y,0:bin_size:0.2);
BINS = EDGES(2:end) -(EDGES(2) - EDGES(1))/2;
h = bar(sfax3_hist,BINS,N,'FaceColor',colors2(3,:)*0.5,'EdgeColor','none','FaceAlpha',0.3);
Xq = BINS;
Vq = movmean(N,5);
hold(sfax3_hist,'on');
plot(sfax3_hist,Xq,Vq,'-r','Color',colors2(3,:),'LineWidth',1);
set(sfax3_hist,'box','off','YColor','w','box','on');
ylim(sfax3_hist,[0,530])
xlabel(sfax3_hist,'L2-norm');

sfax4 = axes(sfig);
sfax4.Units = 'normalized';
sfax4.Position = [0.78,0.65,0.15,0.3];
annotation(sfig, 'textbox', [0.75, 0.99, 0.01, 0.01], 'string', 'd', 'edgecolor', 'none','fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
temp_map = bin_maps.EDR';
imagesc(sfax4,temp_map);colormap(sfax4,cmap);axis(sfax4,'image','xy');
hcb = colorbar(sfax4);
hcb.Label.String = sprintf('%s','Bin counts');
sfax4.Position = [0.78,0.65,0.15,0.3];
set(sfax4,'XTick',1,'XTickLabel','')
set(sfax4,'YTick',(11:10:length(ybin))-0.5,'YTickLabel',ybin(11:10:length(ybin)))
xlabel('L2-norm');
set(sfax4,'FontSize',10,'FontWeight','n')
title('EDR','FontSize',10,'FontWeight','n')
cr = corr(l2_err.EDR(:,1), abs(l2_err.EDR(:,2)), 'Rows', 'complete');
text(sfax4,2,48,sprintf('r = %.03f',cr),'FontSize',10,'FontWeight','n');

sfax4_brain = axes(sfig);
sfax4_brain.Units = 'normalized';
sfax4_brain.Position = [0.85,0.81,0.06,0.12];
cdata = eigenmodes_l2norm.EDR;
patch(sfax4_brain, 'Vertices', surface_midthickness.vertices, 'Faces', surface_midthickness.faces, 'FaceVertexCData', cdata, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap(sfax4_brain,cmap);
caxis(sfax4_brain,[0.02,0.15]);
material(sfax4_brain,'dull');
axis(sfax4_brain,'vis3d');
axis(sfax4_brain,'off');
axis(sfax4_brain,'image');
view(sfax4_brain,270,0);
camlight(sfax4_brain,'headlight')

sfax4_hist = axes(sfig);
sfax4_hist.Units = 'normalized';
sfax4_hist.Position = [0.78,0.58,0.15,0.06];
Y = eigenmodes_l2norm.EDR;
[N,EDGES] = histcounts(Y,0:bin_size:0.2);
BINS = EDGES(2:end) -(EDGES(2) - EDGES(1))/2;
h = bar(sfax4_hist,BINS,N,'FaceColor',colors2(4,:)*0.5,'EdgeColor','none','FaceAlpha',0.3);
Xq = BINS;
Vq = movmean(N,5);
hold(sfax4_hist,'on');
plot(sfax4_hist,Xq,Vq,'-r','Color',colors2(4,:),'LineWidth',1);
set(sfax4_hist,'box','off','YColor','w','box','on');
ylim(sfax4_hist,[0,450])
xlabel(sfax4_hist,'L2-norm');

% Save
set(sfig,'Resize','on','PaperPositionMode','auto','PaperUnits','points','PaperSize',sfig.Position([3,4]) + 1);drawnow;
% saveas(sfig,'./suppl_figure6.pdf');
