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

cfg = []
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
table_SWandBehav=[];
table_P2PandBehav=[];
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
    %load([data_path filesep 'Preproc' filesep 'fixThr_CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    % 1: Subject Number
    % 2: Block Number
    % 3: Electrode Number
    % 4: P2P amplitude
    % 5: Start slow wave (sample from block onset)
    % 12: Downward Slope
    % 13: Upward Slope
    nFc=nFc+1;
    temp_SWandBehav=[];
    temp_P2PandBehav=[];
    for nBl=1:8
        these_SWelectrodes=slow_Waves(slow_Waves(:,2)==nBl,3);
        nout=hist(these_SWelectrodes,1:64);
        nout=nout(matching_elec);

        sub_table_behav=table_behav(table_behav.BlockN==nBl,:);
        duration_block=(sub_table_behav.Sample(end)-sub_table_behav.Sample(1))/hdr.Fs/60;
        densSW=nout/duration_block;
        all_slowWaves(nFc,nBl,:)=densSW;
        
       temp_SWandBehav=[temp_SWandBehav ; [nFc nBl nanmean(table_behav.corrTG(table_behav.BlockN==nBl & table_behav.StimType==1))  nanmean(table_behav.corrNT(table_behav.BlockN==nBl & table_behav.StimType==0))  nanmean(table_behav.RT(table_behav.BlockN==nBl  & table_behav.StimType==1)) ...
           densSW nanmean(densSW)]];
        these_Amplitude=slow_Waves(slow_Waves(:,2)==nBl,4);
        byElec_Amplitude=nan(1,64);
        for nE=1:64
            byElec_Amplitude(nE)=mean(these_Amplitude(these_SWelectrodes==nE));
        end
         all_slowWaves_P2P(nFc,nBl,:)=byElec_Amplitude(matching_elec);
         
         temp_P2PandBehav=[temp_P2PandBehav ; [nFc nBl nanmean(table_behav.corrTG(table_behav.BlockN==nBl & table_behav.StimType==1))  nanmean(table_behav.corrNT(table_behav.BlockN==nBl & table_behav.StimType==0))  nanmean(table_behav.RT(table_behav.BlockN==nBl  & table_behav.StimType==1)) ...
           byElec_Amplitude(matching_elec) nanmean(byElec_Amplitude(matching_elec))]];
       
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
    
    table_SWandBehav=[table_SWandBehav ; mean(temp_SWandBehav,1)];
    table_P2PandBehav=[table_P2PandBehav ; mean(temp_P2PandBehav,1)];
       
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
% figure;
% hp=[];
% [~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves(:,:,match_str(ChanLabels,'Oz'))),0,'k',0,'-',0.5,1,0,1,2);
% title('Averaged slow-waves density across blocks at Oz');
% 
% % ADHD and controls separately
% 
% figure;
% hp=[];
% [~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz'))),0,'b',0,'-',0.5,1,0,1,2);
% hold on;
% [~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz'))),0,'r',0,'-',0.5,1,0,1,2);
% hold on;
% legend(hp,{'Controls','ADHDs'})
% title('Averaged slow-waves density across blocks at Fz');

%RaincloudPlots density
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
h1 = raincloud_plot(squeeze(nanmean(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')),2)), 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',2,'bound_data',[0 100]);
h2 = raincloud_plot(squeeze(nanmean(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')),2)), 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',2,'bound_data',[0 100]);
set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
%set(gca,'XLim', [-30 40], 'YLim', ylim.*[0 0.05]);
format_fig; title('Slow Waves density for Controls and ADHDs'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

[h, pV_diffGroup,~,stat_diffGroup]=ttest2(squeeze(mean(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')),2)),squeeze(mean(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')),2)));
fprintf('... unpaired t-test between groups on SW density on Fz : p=%g, t-value=%g, df=%g\n',...
    pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df) % t(18)=0.90, p=0.38

%% Scalp topographies of average SW density

% all subjects
SW_topo=squeeze(nanmean(nanmean(all_slowWaves(:,:,ismember(ChanLabels,layout.label)),1),2));
figure;
%subplot(1,3,1);
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW density')
caxis([5 11.5])

% ADHD and controls separately
%subplot(1,3,2);
figure;
SW_topo=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar
title('Topography of average SW density for Controls')
caxis([5 11.5])

figure;
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
for nE=1:size(all_slowWaves,3)
    A=squeeze(nanmean(all_slowWaves(match_str(group_SW,'Control'),:,nE),2));
    B=squeeze(nanmean(all_slowWaves(match_str(group_SW,'ADHD'),:,nE),2));
   [h,pV,~,stat]=ttest2(B,A); 
   temp_topo_tval(nE)=stat.tstat;
   temp_topo_pval(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference SW density ADHD/Control (tvalue)')
caxis([-1 1]*4)
if ~isempty(find(temp_topo_pval<0.05))
ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

%% 
% Average P2P amplitude across blocks at Cz
% all subjects
% figure;
% hp=[];
% [~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(:,:,match_str(ChanLabels,'Fz'))),0,'k',0,'-',0.5,1,0,1,2);
% title('Averaged slow-waves P2P across blocks at Fz');
% 
% % ADHD and controls separately
% 
% figure;
% hp=[];
% [~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz'))),0,'b',0,'-',0.5,1,0,1,2);
% hold on;
% [~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz'))),0,'r',0,'-',0.5,1,0,1,2);
% hold on;
% legend(hp,{'Controls','ADHDs'})
% title('Averaged slow-waves P2P across blocks at Fz');

%RaincloudPlots P2P amplitude
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
h1 = raincloud_plot(squeeze(mean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')),2)), 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',1,'bound_data',[0 100]);
h2 = raincloud_plot(squeeze(mean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')),2)), 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',1,'bound_data',[0 100]);
set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
%set(gca,'XLim', [-10 50], 'YLim', ylim.*[0 0.2]);
format_fig; title('Slow Waves P2P amplitude for Controls and ADHDs'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

[h, pV_diffGroup,~,stat_diffGroup]=ttest2(squeeze(mean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Fz')),2)),squeeze(mean(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Fz')),2)));
fprintf('... unpaired t-test between groups on SW P2P amplitude on Fz : p=%g, t-value=%g, df=%g\n',...
    pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df) % t(18)=0.90, p=0.38

% Scalp topographies of P2P amplitude
% all subjects

figure;
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(:,:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW P2P')
caxis([20 50])

% ADHD and controls separately
figure;
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar
title('Topography of average SW P2P for Controls')
caxis([20 50])

figure;
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,ismember(ChanLabels,layout.label)))));
simpleTopoPlot_ft(SW_topo, layout,'on',[],0,1);
colorbar;
title('Topography of average SW P2P for ADHDs')
caxis([20 50])

%T-values P2P amplitude topo
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval_P2P=[];
temp_topo_pval_P2P=[];
for nE=1:size(all_slowWaves_P2P,3)
    A=squeeze(nanmean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,nE),2));
    B=squeeze(nanmean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,nE),2));
   [h,pV,~,stat]=ttest2(B,A); 
   temp_topo_tval_P2P(nE)=stat.tstat;
   temp_topo_pval_P2P(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval_P2P,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference SW amplitude ADHD/Control (tvalue)')
caxis([-1 1]*4)
if ~isempty(find(temp_topo_pval_P2P<0.05))
ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_P2P<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end
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

%T-VALUE SW SLOPE
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval_DWS=[];
temp_topo_pval_DWS=[];
for nE=1:size(all_slowWaves_DWslope,3)
    A=squeeze(nanmean(all_slowWaves_DWslope(match_str(group_SW,'Control'),:,nE),2));
    B=squeeze(nanmean(all_slowWaves_DWslope(match_str(group_SW,'ADHD'),:,nE),2));
   [h,pV,~,stat]=ttest2(B,A); 
   temp_topo_tval_DWS(nE)=stat.tstat;
   temp_topo_pval_DWS(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval_DWS,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference DW slope ADHD/Control (tvalue)')
caxis([-1 1]*4)
if ~isempty(find(temp_topo_pval_DWS<0.05))
ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_DWS<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end
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

%T-VALUE SW SLOPE
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval_UWS=[];
temp_topo_pval_UWS=[];
for nE=1:size(all_slowWaves_UPWslope,3)
    A=squeeze(nanmean(all_slowWaves_UPWslope(match_str(group_SW,'Control'),:,nE),2));
    B=squeeze(nanmean(all_slowWaves_UPWslope(match_str(group_SW,'ADHD'),:,nE),2));
   [h,pV,~,stat]=ttest2(B,A); 
   temp_topo_tval_UWS(nE)=stat.tstat;
   temp_topo_pval_UWS(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval_UWS,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference UW slope ADHD/Control (tvalue)')
caxis([-1 1]*4)
if ~isempty(find(temp_topo_pvalUWS<0.05))
ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_UWS<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

%%
% Repeated Measures ANOVA
t = table(group_SW,all_slowWaves(:,1,:),all_slowWaves(:,2,:),all_slowWaves(:,3,:),all_slowWaves(:,4,:),all_slowWaves(:,5,:),all_slowWaves(:,6,:),all_slowWaves(:,7,:),all_slowWaves(:,8,:),...
'VariableNames',{'Group','B1','B2','B3','B4','B5','B6','B7','B8'});
Time = [1 2 3 4 5 6 7 8]';
rm = fitrm(t,'B1-B8 ~ Group','WithinDesign',Time,'WithinModel','orthogonalcontrasts')
ranovatbl = ranova(rm)

%% P2P Amplitude and behav
Colors2=[247,104,161;
    212,185,218;
    251,128,114;
    255,255,179]/256;

figure;
subplot(1,3,1);
simpleCorPlot(table_P2PandBehav(:,end),1-table_P2PandBehav(:,3),{'o',Colors2(2,:),Colors2(2,:),50},'Spearman');
title('Correlation Miss/SW P2P');


subplot(1,3,2);
simpleCorPlot(table_P2PandBehav(:,end),1-table_P2PandBehav(:,4),{'o',Colors2(1,:),Colors2(1,:),50},'Spearman');
title('Correlation FA/SW P2P')

subplot(1,3,3);
simpleCorPlot(table_P2PandBehav(:,end),table_P2PandBehav(:,5),{'o',Colors2(3,:),Colors2(3,:),50},'Spearman');
title('Correlation RT/SW P2P')

myFigPos=[182         369        1150         428]
set(gcf,'Position',myFigPos);

%% Topo correlation P2P amplitude and behav
%Miss
temp_topo_rval_P2PMiss=[];
temp_topo_pval_P2PMiss=[];
for nE=1:size(all_slowWaves_P2P,3)
    [r pV]=corr(table_P2PandBehav(:,5+nE),1-table_P2PandBehav(:,3),'rows','pairwise','type','Spearman');
    temp_topo_rval_P2PMiss(nE)=r;
    temp_topo_pval_P2PMiss(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_rval_P2PMiss,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography Correlation coeff Miss/SW P2P')
caxis([-1 1]*1)
if ~isempty(find(temp_topo_pval_P2PMiss<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_P2PMiss<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

%FA
temp_topo_rval_P2PFA=[];
temp_topo_pval_P2PFA=[];
for nE=1:size(all_slowWaves_P2P,3)
    [r pV]=corr(table_P2PandBehav(:,5+nE),1-table_P2PandBehav(:,4),'rows','pairwise','type','Spearman');
    temp_topo_rval_P2PFA(nE)=r;
    temp_topo_pval_P2PFA(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_rval_P2PFA,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography Correlation coeff FA/SW P2P')
caxis([-1 1]*1)
if ~isempty(find(temp_topo_pval_P2PFA<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_P2PFA<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

%RT
temp_topo_rval_P2PRT=[];
temp_topo_pval_P2PRT=[];
for nE=1:size(all_slowWaves_P2P,3)
    [r pV]=corr(table_P2PandBehav(:,5+nE),table_P2PandBehav(:,5),'rows','pairwise','type','Spearman');
    temp_topo_rval_P2PRT(nE)=r;
    temp_topo_pval_P2PRT(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_rval_P2PRT,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography Correlation coeff RT/SW P2P')
caxis([-1 1]*1)
if ~isempty(find(temp_topo_pval_P2PRT<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_P2PRT<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end
%% SW density and behav
figure;
subplot(1,3,1);
simpleCorPlot(table_SWandBehav(:,end),1-table_SWandBehav(:,3),{'o',Colors2(2,:),Colors2(2,:),50},'Spearman');
title('Correlation Miss/SW density');


subplot(1,3,2);
simpleCorPlot(table_SWandBehav(:,end),1-table_SWandBehav(:,4),{'o',Colors2(1,:),Colors2(1,:),50},'Spearman');
title('Correlation FA/SW density')

subplot(1,3,3);
simpleCorPlot(table_SWandBehav(:,end),table_SWandBehav(:,5),{'o',Colors2(3,:),Colors2(3,:),50},'Spearman');
title('Correlation RT/SW density')

myFigPos=[182         369        1150         428]
set(gcf,'Position',myFigPos);
%% Topo correlation SW density and behav
%Miss
temp_topo_rval_DMiss=[];
temp_topo_pval_DMiss=[];
for nE=1:size(all_slowWaves,3)
    [r pV]=corr(table_SWandBehav(:,5+nE),1-table_SWandBehav(:,3),'rows','pairwise','type','Spearman');
    temp_topo_rval_DMiss(nE)=r;
    temp_topo_pval_DMiss(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_rval_DMiss,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography Correlation coeff Miss/SW Density')
caxis([-1 1]*1)
if ~isempty(find(temp_topo_pval_DMiss<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_DMiss<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

%FA
temp_topo_rval_DFA=[];
temp_topo_pval_DFA=[];
for nE=1:size(all_slowWaves,3)
    [r pV]=corr(table_SWandBehav(:,5+nE),1-table_SWandBehav(:,4),'rows','pairwise','type','Spearman');
    temp_topo_rval_DFA(nE)=r;
    temp_topo_pval_DFA(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_rval_DFA,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography Correlation coeff FA/SW Density')
caxis([-1 1]*1)
if ~isempty(find(temp_topo_pval_DFA<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_DFA<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

%RT
temp_topo_rval_DRT=[];
temp_topo_pval_DRT=[];
for nE=1:size(all_slowWaves,3)
    [r pV]=corr(table_SWandBehav(:,5+nE),table_SWandBehav(:,5),'rows','pairwise','type','Spearman');
    temp_topo_rval_DRT(nE)=r;
    temp_topo_pval_DRT(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_rval_DRT,layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography Correlation coeff RT/SW Density')
caxis([-1 1]*1)
if ~isempty(find(temp_topo_pval_DRT<0.05))
    ft_plot_lay_me(layout, 'chanindx',find(temp_topo_pval_DRT<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end