%%
clear all
close all

run ../localdef_ADHD_CTET.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep '*' filesep '*' filesep '*CTET*.bdf']);

cfg = [];
cfg.layout = 'biosemi64.lay';
layout=ft_prepare_layout(cfg);

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel = layout.label;
cfg.channel(match_str(layout.label,{'Iz','P7','P8'}))=[];
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

ChanLabels=layout.label(1:end-2);
%%
all_slowWaves=[];
nFc=0;
nFc4=0;
group_SW=[];
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    if exist([save_path filesep 'CTET_ADHD_behav_' file_name(1:end-4) '.txt'])==0
        warning(sprintf('missing behavioural file for %s\n',file_name(1:end-4)));
        continue;
    end
    table_behav=readtable([save_path filesep 'CTET_ADHD_behav_' file_name(1:end-4) '.txt']);
    hdr=ft_read_header([folder_name filesep file_name]);
    matching_elec=[];
    for nE=1:length(layout.label)-2
        matching_elec(nE)=(match_str(hdr.label,layout.label(nE)));
    end
        
    %load([data_path filesep 'Preproc' filesep 'CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    load([data_path filesep 'Preproc' filesep 'relThrCTR_CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    % 1: Subject Number
    % 2: Block Number
    % 3: Electrode Number
    % 4: P2P amplitude
    % 5: Start slow wave (sample from block onset)
    % 12: Downward Slope
    % 13: Upward Slope
    nFc=nFc+1;
    for nBl=1:8
        these_SWelectrodes=slow_Waves(slow_Waves(:,2)==nBl,3);
        nout=hist(these_SWelectrodes,1:64);
        nout=nout(matching_elec);

        sub_table_behav=table_behav(table_behav.BlockN==nBl,:);
        duration_block=(sub_table_behav.Sample(end)-sub_table_behav.Sample(1))/hdr.Fs/60;
        densSW=nout/duration_block;
        all_slowWaves(nFc,nBl,:)=densSW;
        
        these_Amplitude=slow_Waves(slow_Waves(:,2)==nBl,4);
        byElec_Amplitude=nan(1,64);
        for nE=1:64
            byElec_Amplitude(nE)=mean(these_Amplitude(these_SWelectrodes==nE));
        end
         all_slowWaves_P2P(nFc,nBl,:)=byElec_Amplitude(matching_elec);
         
        these_DWslope=slow_Waves(slow_Waves(:,2)==nBl,12);
        byElec_Amplitude=nan(1,64);
        for nE=1:64
            byElec_Amplitude(nE)=mean(these_DWslope(these_SWelectrodes==nE));
        end
         all_slowWaves_DWslope(nFc,nBl,:)=byElec_Amplitude(matching_elec);
         
        these_UPWslope=slow_Waves(slow_Waves(:,2)==nBl,13);
        byElec_Amplitude=nan(1,64);
        for nE=1:64
            byElec_Amplitude(nE)=mean(these_UPWslope(these_SWelectrodes==nE));
        end
         all_slowWaves_UPWslope(nFc,nBl,:)=byElec_Amplitude(matching_elec);
   end
    
    orifoldername=files(nF).folder;
    if isempty(findstr(orifoldername,'controls'))==0
        group_SW{nFc}='Control';
    elseif isempty(findstr(orifoldername,'adhds'))==0
        nFc4=nFc4+1;
        group_SW{nFc}='ADHD';
    end
end

%%
% Average SW density across blocks
% all subjects
figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves(:,:,match_str(ChanLabels,'Oz'))),0,'k',0,'-',0.5,1,0,1,2);
title('Averaged slow-waves density across blocks at Oz');

% ADHD and controls separately

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz'))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('Averaged slow-waves density across blocks at Fz');

%RaincloudPlots density
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
h1 = raincloud_plot(squeeze(mean(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')))), 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.4,'bound_data',[0 100]);
h2 = raincloud_plot(squeeze(mean(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')))), 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.4,'bound_data',[0 100]);
set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
%set(gca,'XLim', [-30 40], 'YLim', ylim.*[0 0.05]);
format_fig; title('Slow Waves density for Controls and ADHDs'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

[h, pV_diffGroup,~,stat_diffGroup]=ttest2(squeeze(mean(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')))),squeeze(mean(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')))));
fprintf('... unpaired t-test between groups on SW density on Fz : p=%g, t-value=%g, df=%g\n',...
    pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df) % t(18)=0.90, p=0.38

%% Scalp topographies of average SW density

% all subjects
SW_topo=squeeze(nanmean(nanmean(all_slowWaves(:,:,ismember(ChanLabels,layout.label)),1),2));
figure;
subplot(1,3,1);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW density')
caxis([5 11.5])

% ADHD and controls separately
subplot(1,3,2);
SW_topo=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar
title('Topography of average SW density for Controls')
caxis([5 11.5])

subplot(1,3,3);
SW_topo=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW density for ADHDs')
caxis([5 11.5])

%T-values density topo
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(all_slowWaves,2)
    A=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
    B=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
   [h,pV,~,stat]=ttest2(A,B); 
   temp_topo_tval(nE)=stat.tstat;
   temp_topo_pval(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval(matching_elec),layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference SW density ADHD/Control (tvalue)')
%caxis([-1 1]*4)

%% 
% Average P2P amplitude across blocks at Cz
% all subjects
figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(:,:,match_str(ChanLabels,'Fz'))),0,'k',0,'-',0.5,1,0,1,2);
title('Averaged slow-waves P2P across blocks at Fz');

% ADHD and controls separately

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz'))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('Averaged slow-waves P2P across blocks at Fz');

%RaincloudPlots P2P amplitude
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
h1 = raincloud_plot(squeeze(mean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')))), 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.6,'bound_data',[0 100]);
h2 = raincloud_plot(squeeze(mean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')))), 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.6,'bound_data',[0 100]);
set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
%set(gca,'XLim', [-10 50], 'YLim', ylim.*[0 0.2]);
format_fig; title('Slow Waves P2P amplitude for Controls and ADHDs'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

[h, pV_diffGroup,~,stat_diffGroup]=ttest2(squeeze(mean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')))),squeeze(mean(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')))));
fprintf('... unpaired t-test between groups on SW P2P amplitude on Fz : p=%g, t-value=%g, df=%g\n',...
    pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df) % t(18)=0.90, p=0.38

% Scalp topographies of P2P amplitude
% all subjects

figure;
subplot(1,3,1);
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(:,:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW P2P')
caxis([20 50])

% ADHD and controls separately
subplot(1,3,2);
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar
title('Topography of average SW P2P for Controls')
caxis([20 50])

subplot(1,3,3);
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW P2P for ADHDs')
caxis([20 50])

%T-values P2P amplitude topo
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(all_slowWaves_P2P,2)
    A=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
    B=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
   [h,pV,~,stat]=ttest2(A,B); 
   temp_topo_tval(nE)=stat.tstat;
   temp_topo_pval(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval(matching_elec),layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference SW P2P amplitude ADHD/Control (tvalue)')
%caxis([-1 1]*4)
%%
% Average Downward slope across blocks
% all subjects
figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_DWslope(:,:,match_str(ChanLabels,'Oz'))),0,'k',0,'-',0.5,1,0,1,2);
title('Averaged slow-waves Downward slope across blocks at Oz');

% ADHD and controls separately

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_DWslope(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Oz'))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves_DWslope(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Oz'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('Averaged slow-waves Downward slope across blocks at Oz');

% Scalp topographies of DW slope
% all subjects

SW_topo=squeeze(nanmean(nanmean(all_slowWaves_DWslope(:,:,ismember(ChanLabels,layout.label)))));
subplot(1,3,1);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW Downward slope')
caxis([190 500])

% ADHD and controls separately
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_DWslope(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
subplot(1,3,2);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar
title('Topography of average SW Downward slope for Controls')
caxis([190 500])

SW_topo=squeeze(nanmean(nanmean(all_slowWaves_DWslope(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
subplot(1,3,3);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW Downward slope for ADHDs')
caxis([190 500])

%%
% Average Upward slope across blocks
% all subjects
figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_UPWslope(:,:,match_str(ChanLabels,'Oz'))),0,'k',0,'-',0.5,1,0,1,2);
title('Averaged slow-waves Upward slope across blocks at Oz');

% ADHD and controls separately

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_UPWslope(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Oz'))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves_UPWslope(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Oz'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('Averaged slow-waves Upward slope across blocks at Oz');

% Scalp topographies of DW slope
% all subjects
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_UPWslope(:,:,ismember(ChanLabels,layout.label)))));
figure;
subplot(1,3,1);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW Upward slope')
caxis([190 500])

% ADHD and controls separately
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_UPWslope(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
subplot(1,3,2);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar
title('Topography of average SW Upward slope for Controls')
caxis([190 500])

SW_topo=squeeze(nanmean(nanmean(all_slowWaves_UPWslope(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
subplot(1,3,3);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW Upward slope for ADHDs')
caxis([190 500])

%%
% Repeated Measures ANOVA
t = table(group_SW,all_slowWaves(:,1,:),all_slowWaves(:,2,:),all_slowWaves(:,3,:),all_slowWaves(:,4,:),all_slowWaves(:,5,:),all_slowWaves(:,6,:),all_slowWaves(:,7,:),all_slowWaves(:,8,:),...
'VariableNames',{'Group','B1','B2','B3','B4','B5','B6','B7','B8'});
Time = [1 2 3 4 5 6 7 8]';
rm = fitrm(t,'B1-B8 ~ Group','WithinDesign',Time,'WithinModel','orthogonalcontrasts')
ranovatbl = ranova(rm)