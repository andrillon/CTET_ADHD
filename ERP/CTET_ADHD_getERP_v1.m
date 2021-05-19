%%
clear all;
% close all;
run ../localdef_ADHD_CTET.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep 'Preproc' filesep 'CIcfe_ft_*.mat']);

% table=readtable([save_path 'CTET_behav_res.txt']);
% load('../ICA_Artifacts.mat')

%% Layout
load('../ADHD_CTET_task_BadChannels.mat')
redo=0;
%%
nFc=0;
nFc4=0;
all_ERP_NT=[];
all_ERP_TG=[];

all_ERP_NT_offset=[];
all_ERP_TG_offset=[];
group_PowDataEO=[];
    
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubIDlong=file_name(10:end-4);
    SubID=file_name(10:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(end)-1);
    if exist([save_path filesep 'CTET_ADHD_behav_' SubID '_CTET.txt'])==0
        continue;
    end
    table=readtable([save_path filesep 'CTET_ADHD_behav_' SubID '_CTET.txt']);
    orifile=dir([data_path filesep '*' filesep  '*' filesep SubIDlong '.bdf']);

 
    nFc=nFc+1;
    fprintf('... processing %s\n',file_name);
    load([data_path filesep  'Preproc'  filesep file_name(1:end-4)]);
    
    cfg=[];
    cfg.reref           = 'yes';
    cfg.refchannel      = 'all';
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.2 0];
    data_clean = ft_preprocessing(cfg,data);
    
    %%% take out trial
    thisF=match_str(badChannels_table.FileName,['fe_ft_' SubIDlong]);
    if isempty(thisF)
        continue;
    end
    eval(['badTrials=[' badChannels_table.Bad_Trials{thisF} '];']);
    table(ismember(table.TrialN,badTrials),:)=[];
    table.StimType(find(diff(table.BlockN)==1 | diff(table.StimType)==1))=NaN;
    
    % COMPUTE BOTH ONSET AND OFFSET ERP
    cfgerp        = [];
    %     cfgerp.latency        = [-0.2 0.8];
    cfgerp.trials = find(table.StimType==0); %cfgerp.trials(cfgerp.trials>length(data_clean.trial))=[];
    av_data_NT = ft_timelockanalysis(cfgerp, data_clean);
    cfgerp.trials = find(table.StimType==1); %cfgerp.trials(cfgerp.trials>length(data_clean.trial))=[];
    av_data_TG = ft_timelockanalysis(cfgerp, data_clean);
    
    ERP_NT=av_data_NT.avg; %(:,av_data_NT.time>-0.2 & av_data_NT.time<0.8);
    all_ERP_NT(nFc,:,:)=ERP_NT;
    
    ERP_TG=av_data_TG.avg; %(:,av_data_TG.time>-0.2 & av_data_TG.time<0.8);
    all_ERP_TG(nFc,:,:)=ERP_TG;
        
    cfgerp2        = [];
    cfgerp2.trials = find(table.StimType==0)+1; cfgerp2.trials(cfgerp2.trials>length(data_clean.trial))=[];
    av_data_NT_offset = ft_timelockanalysis(cfgerp2, data_clean);
    cfgerp2.trials = find(table.StimType==1)+1; cfgerp2.trials(cfgerp2.trials>length(data_clean.trial))=[];
    av_data_TG_offset = ft_timelockanalysis(cfgerp2, data_clean);
    
    ERP_NT_offset=av_data_NT_offset.avg(:,av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0.8)-repmat(mean(av_data_NT_offset.avg(:,av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0),2),1,sum(av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0.8));
    all_ERP_NT_offset(nFc,:,:)=ERP_NT_offset;
    
    ERP_TG_offset=av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8)-repmat(mean(av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0),2),1,sum(av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8));
    all_ERP_TG_offset(nFc,:,:)=ERP_TG_offset;
    
    
    orifoldername=orifile.folder;
    if isempty(findstr(orifoldername,'controls'))==0
        group_PowDataEO{nFc}='Control';
    elseif isempty(findstr(orifoldername,'adhds'))==0
        nFc4=nFc4+1;
        group_PowDataEO{nFc}='ADHD';
    end
end

xTime=av_data_NT.time;
chLabels=av_data_NT.label;
xTime_offset=av_data_TG_offset.time(av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8);

%%
diff_all_ERP=all_ERP_TG-all_ERP_NT;
diff_all_ERP_offset=all_ERP_TG_offset-all_ERP_NT_offset;
%%
%Plots
%Event-related potentials non-target vs target 

thisCh=match_str(chLabels,'POz');

Colors1=[233,163,201;
247,247,247;
161,215,106]/256;

figure;
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NT(:,thisCh,:)),0,Colors1(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_TG(:,thisCh,:)),0,Colors1(3,:),0,'-',0.5,1,0,1,2);
hold on;
% legend(hp,{'Non Target','Target'})
title('Event-related potentials non-target vs target');
xlim([-0.2 1.8])
format_fig;


figure;
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_TG(:,thisCh,:))-squeeze(all_ERP_NT(:,thisCh,:)),0,'r',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);

%%
thisCh=match_str(chLabels,'POz');
diffTG_NT=squeeze(mean(all_ERP_TG(:,thisCh,xTime>1 & xTime<1.8)-all_ERP_NT(:,thisCh,xTime>1 & xTime<1.8),3));

%DoubleDiff = diff_all_ERP(match_str(group_PowDataEO,'Control'))-diff_all_ERP(match_str(group_PowDataEO,'ADHD'))
%squeeze(mean(mean(diffTG_NT((match_str(group_PowDataEO,'Control')),thisCh,xTime>1.1 & xTime<1.8)-diffTG_NT((match_str(group_PowDataEO,'ADHD'),thisCh,xTime>1.1 & xTime<1.8),3));
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
h1 = raincloud_plot(diffTG_NT(match_str(group_PowDataEO,'Control')), 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',8);
h2 = raincloud_plot(diffTG_NT(match_str(group_PowDataEO,'ADHD')), 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',8);

[h, pV_diffGroup,~,stat_diffGroup]=ttest2(diffTG_NT(match_str(group_PowDataEO,'Control')),diffTG_NT(match_str(group_PowDataEO,'ADHD')));
fprintf('... unpaired t-test between groups on diff ERP between TG and NT on Pz for 1.1-1.8s: p=%g, t-value=%g, df=%g\n',...
    pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df) % t(18)=0.90, p=0.38

[h, pV_CTR,~,stat_CTR]=ttest(diffTG_NT(match_str(group_PowDataEO,'Control')),0);
fprintf('... one-sample t-test with 0 for Controls on diff ERP between TG and NT on Pz for 1.1-1.8s: p=%g, t-value=%g, df=%g\n',...
    pV_CTR,stat_CTR.tstat,stat_CTR.df)
[h, pV_ADHD,~,stat_ADHD]=ttest(diffTG_NT(match_str(group_PowDataEO,'ADHD')),0);
fprintf('... unpaired t-test with 0 for ADHD on diff ERP between TG and NT on Pz for 1.1-1.8s: p=%g, t-value=%g, df=%g\n',...
    pV_ADHD,stat_ADHD.tstat,stat_ADHD.df)

% %%
% %Event-related potentials non-target vs target OFFSET
% thisCh=match_str(chLabels,'Fz');
% figure;
% hp=[];
% [~,hp(1)]=simpleTplot(xTime_offset,squeeze(all_ERP_NT_offset(:,thisCh,:)),0,'b',0,'-',0.5,1,0,1,2);
% hold on;
% [~,hp(2)]=simpleTplot(xTime_offset,squeeze(all_ERP_TG_offset(:,thisCh,:)),0,'r',0,'-',0.5,1,0,1,2);
% hold on;
% legend(hp,{'Non Target','Target'})
% title('Event-related potentials non-target vs target OFFSET');
% format_fig;

%%
%%% Plot topography [0.1-0.3]s post-offset
cfg = [];
cfg.layout = 'biosemi64.lay';
layout=ft_prepare_layout(cfg);

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel = layout.label;
cfg.channel(match_str(layout.label,{'Iz','P7','P8'}))=[];
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

matching_elec=[];
    for nE=1:length(layout.label)-2
        matching_elec(nE)=(match_str(data.label,layout.label(nE)));
    end
    
% %Non-target trials offset-locked
% temp_topo=squeeze(mean(mean(all_ERP_NT_offset(:,:,xTime_offset>0.1 & xTime_offset<0.3),3),1));
% figure;
% simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
% colorbar;
% title('Topography non-target trials [0.1-0.3]s post-offset')
% 
% %Target trials offset-locked
% temp_topo=squeeze(mean(mean(all_ERP_TG_offset(:,:,xTime_offset>0.1 & xTime_offset<0.3),3),1));
% figure;
% simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
% colorbar;
% title('Topography target trials [0.1-0.3]s post-offset')

%Difference TG/NT trials onset-locked
temp_topo=squeeze(mean(mean(diff_all_ERP(:,matching_elec,xTime>1.1 & xTime<1.8),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colorbar;
title('Topography difference target/non target trials [1.1-1.8]s post-onset')
caxis([-1 1]*5.2)

temp_topo_CTR=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),matching_elec,xTime>1.1 & xTime<1.8),3),1));
figure;
simpleTopoPlot_ft(temp_topo_CTR', layout,'on',[],0,1);
colorbar;
title('Topography difference target/non target trials [1.1-1.8]s post-onset for Controls')
caxis([-1 1]*5.2)

temp_topo_ADHD=squeeze(mean(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime>1.1 & xTime<1.8),3),1));
figure;
simpleTopoPlot_ft(temp_topo_ADHD', layout,'on',[],0,1);
colorbar;
title('Topography difference target/non target trials [1.1-1.8]s post-onset for ADHDs')
caxis([-1 1]*5.2)

cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(diff_all_ERP,2)
    A=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'Control'),nE,xTime>1.1 & xTime<1.8),3));
    B=squeeze(mean(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),nE,xTime>1.1 & xTime<1.8),3));
   [h,pV,~,stat]=ttest2(A,B); 
   temp_topo_tval(nE)=stat.tstat;
   temp_topo_pval(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval(matching_elec), layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference target/non target trials [1.1-1.8]s post-onset (tvalue)')
caxis([-1 1]*4)
%%
%%% Plot topography [0.3-0.8]s post-offset
%Non-target trials offset-locked
temp_topo=squeeze(mean(mean(all_ERP_NT_offset(:,matching_elec,xTime_offset>0.3 & xTime_offset<0.8),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography non-target trials [0.3-0.8]s post-offset')

%Target trials offset-locked
temp_topo=squeeze(mean(mean(all_ERP_TG_offset(:,matching_elec,xTime_offset>0.3 & xTime_offset<0.8),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography target trials [0.3-0.8]s post-offset')

%Difference TG/NT trials offset-locked
temp_topo=squeeze(mean(mean(diff_all_ERP_offset(:,matching_elec,xTime_offset>0.3 & xTime_offset<0.8),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography difference target/non target trials [0.3-0.8]s post-offset')

%%
%%Difference TG/NG for all subjects

thisCh=match_str(chLabels,'Oz');
figure;
subplot(1,2,1);
[~,~]=simpleTplot(xTime,squeeze(diff_all_ERP(:,thisCh,:)),0,'k',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
title('ERP differences between Target and Non-target trials for all subjects ONSET');
subplot(1,2,2);
[~,~]=simpleTplot(xTime_offset,squeeze(diff_all_ERP_offset(:,thisCh,:)),0,'k',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
title('ERP differences between Target and Non-target trials for all subjects OFFSET');

%%
%Difference TG/NG for Controls and ADHDs
thisCh=match_str(chLabels,'Fz');

figure;
subplot(1,2,1);
hp=[];
% [~,hp(1)]=simpleTplot(xTime,squeeze(diff_all_ERP(:,thisCh,:)),0,'k');
% hold on;
[~,hp(1)]=simpleTplot(xTime,squeeze(diff_all_ERP(match_str(group_PowDataEO,'Control'),thisCh,:)),0,'b',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(diff_all_ERP(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,'r',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('ERP differences between Target and Non-target trials ONSET');

subplot(1,2,2);
hp=[];
% [~,hp(1)]=simpleTplot(xTime_offset,squeeze(diff_all_ERP_offset(:,thisCh,:)),0,'k');
% hold on;
[~,hp(1)]=simpleTplot(xTime_offset,squeeze(diff_all_ERP_offset(match_str(group_PowDataEO,'Control'),thisCh,:)),0,'b',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime_offset,squeeze(diff_all_ERP_offset(match_str(group_PowDataEO,'ADHD'),thisCh,:)),0,'r',[2 0.05 0.05 1000],'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('ERP differences between Target and Non-target trials OFFSET');

%%
%Topography Controls and ADHD
%%% Plot topography [0.1-0.3]s post-offset
%Non-target trials offset-locked Controls
temp_topo=squeeze(mean(mean(all_ERP_NT_offset(match_str(group_PowDataEO,'Control'),matching_elec,xTime_offset>0.1 & xTime_offset<0.3),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography Controls non-target trials [0.1-0.3]s post-offset')

%Target trials offset-locked Controls
temp_topo=squeeze(mean(mean(all_ERP_TG_offset(match_str(group_PowDataEO,'Control'),matching_elec,xTime_offset>0.1 & xTime_offset<0.3),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography Controls target trials [0.1-0.3]s post-offset')

%Difference TG/NT trials offset-locked Controls
temp_topo=squeeze(mean(mean(diff_all_ERP_offset(match_str(group_PowDataEO,'Control'),matching_elec,xTime_offset>0.1 & xTime_offset<0.3),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography difference target/non target trials [0.1-0.3]s post-offset for Controls')

%Non-target trials offset-locked ADHDs
temp_topo=squeeze(mean(mean(all_ERP_NT_offset(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime_offset>0.1 & xTime_offset<0.3),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography ADHD non-target trials [0.1-0.3]s post-offset')

%Target trials offset-locked ADHDS
temp_topo=squeeze(mean(mean(all_ERP_TG_offset(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime_offset>0.1 & xTime_offset<0.3),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography ADHD target trials [0.1-0.3]s post-offset')

%Difference TG/NT trials offset-locked ADHDs
temp_topo=squeeze(mean(mean(diff_all_ERP_offset(match_str(group_PowDataEO,'ADHD'),matching_elec,xTime_offset>0.1 & xTime_offset<0.3),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Topography difference target/non target trials [0.1-0.3]s post-offset for ADHDs')
%%
%%% Cluster difference ADHD and controls separaterly for onset-locked and
%%% offset-locked data


%%
All_Conds=double(ismember(group_PowDataEO,'Control'))+1;
[realpos_lin realneg_lin]=get_cluster_permutation_aov(squeeze(diff_all_ERP(:,thisCh,:)),All_Conds',...
        0.05,0.05,500,xTime);
    
