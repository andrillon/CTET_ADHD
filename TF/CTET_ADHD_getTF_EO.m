%%
clear all;
close all;

%%
run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;

files=dir([data_path filesep 'Preproc' filesep 'CIcf_ft_*.mat']);


%% loop on subjects
redo=1;
nFc=0;
nFc1=0;
nFc2=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubIDlong=file_name(9:end-4);
    SubID=file_name(9:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    %%% load 
    load([data_path filesep 'Preproc' filesep file_name(1:end-4)]);
    
    %%% Retrieve ADHD or control
    orifile=dir([data_path filesep '*' filesep  '*' filesep SubIDlong '.bdf']);

    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          =  0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
    cfg.toi         =  'all';                         % analysis 2 to 30 Hz in steps of .2 Hz
    cfg.t_ftimwin    =  ones(length(cfg.foi),1).*10;   % length of time window =6 sec
    cfg.keeptrials   = 'yes';
    TFdata           = ft_freqanalysis(cfg, data);
    
    nFc=nFc+1;
    av_PowDataEO(nFc,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
    if isempty(findstr(file_name,'efore'))==0
        nFc1=nFc1+1;
%         av_PowDataEO_Before(nFc1,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
        cond_PowDataEO{nFc}='Before';
    elseif isempty(findstr(file_name,'fter'))==0
        nFc2=nFc2+1;
%         av_PowDataEO_After(nFc2,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
        cond_PowDataEO{nFc}='After';
    end
    
    orifoldername=orifile.folder;
    if isempty(findstr(file_name,'controls'))==0
        group_PowDataEO{nFc}='Control';
    elseif isempty(findstr(file_name,'adhds'))==0
        group_PowDataEO{nFc}='ADHD';
    end
end

%% 1: Plot average power for Cz Oz and Fz
figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Cz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Oz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Fz'),:),1)))

%% 2: Split between before and after 
figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),match_str(data.label,'Fz'),:),1)),'Color','b')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),match_str(data.label,'Fz'),:),1)),'Color','r')
legend({'Before','After'})
title('Fz')

%% 3: Split between ADHD and controls

figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Control'),match_str(data.label,'Fz'),:),1)),'Color','b')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'ADHD'),match_str(data.label,'Fz'),:),1)),'Color','r')
legend({'Controls','ADHD'})
title('Fz')
%% 4: Split between before and after + ADHD and controls 

figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Control'),match_str(cond_PowDataEO,'Before'),match_str(data.label,'Fz'),:),1)),'Color','b')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Control'),match_str(cond_PowDataEO,'After'),match_str(data.label,'Fz'),:),1)),'Color','r')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'Before'),match_str(data.label,'Fz'),:),1)),'Color','g')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'After'),match_str(data.label,'Fz'),:),1)),'Color','y')
legend({'Controls Before','Controls After','ADHD Before','ADHD After'})
title('Fz')

