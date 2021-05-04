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
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
load('../ADHD_CTET_task_BadChannels.mat')

redo=1;
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
%Plots

thisCh=match_str(chLabels,'Pz');

figure;
hp=[];
[~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NT(:,thisCh,:)),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_TG(:,thisCh,:)),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Non Target','Target'})
title('Event-related potentials non-target vs target');


%%
thisCh=match_str(chLabels,'Oz');
figure;
hp=[];
[~,hp(1)]=simpleTplot(xTime_offset,squeeze(all_ERP_NT_offset(:,thisCh,:)),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(xTime_offset,squeeze(all_ERP_TG_offset(:,thisCh,:)),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Non Target','Target'})
title('Event-related potentials offset non-target vs target');

%%
%Difference TG/NG for all subjects
thisCh=match_str(chLabels,'Pz');

figure;
hp=[];
diff=all_ERP_TG-all_ERP_NT;
[~,hp(1)]=simpleTplot(xTime,squeeze(diff(:,thisCh,:)),0,'b');
hold on;
[~,hp(2)]=simpleTplot(xTime,squeeze(diff(intersect(match_str(chLabels,'Pz'),match_str(group_PowDataEO,'Control'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(xTime,squeeze(diff(intersect(match_str(chLabels,'Pz'),match_str(group_PowDataEO,'ADHD'))),0,'y',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'All subjects','Controls','ADHDs'})
title('ERP differences between Target and Non-target trials');




