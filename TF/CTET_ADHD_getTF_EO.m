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
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(9:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    %%% load 
    load([data_path filesep 'Preproc' filesep file_name(1:end-4)]);
    
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          =  0.5:0.5:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
    cfg.toi         =  'all';                         % analysis 2 to 30 Hz in steps of .2 Hz
    cfg.t_ftimwin    =  ones(length(cfg.foi),1).*10;   % length of time window =6 sec
    cfg.keeptrials   = 'yes';
    TFdata           = ft_freqanalysis(cfg, data);
    
    nFc=nFc+1;
    av_PowDataEO(nFc,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
%     if redo==1 || exist([data_path filesep 'Preproc' filesep 'CIcf_ft_' file_name(1:end-4) '.mat'])==0
%                save([data_path filesep 'Preproc' filesep 'CIcf_ft_' file_name(1:end-4)],'data');
% 
%     end
end

%%%% 1: Plot average power for Cz Oz and Fz


%% 2: Split between before and after 

redo=1;
nFc=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(9:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    if findstr(seps,'Before')||findstr(seps,'before')
        dataBefore = data
    end
    
     if findstr(seps,'After')||findstr(seps,'after')
        dataAfter = data
     end
    
    %%% load 
        load([data_path filesep 'Preproc' filesep file_name(1:end-4)]);
   
            cfg              = [];
            cfg.output       = 'pow';
            cfg.channel      = 'all';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.foi          =  0.5:0.5:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
            cfg.toi          =  'all';                         % analysis 2 to 30 Hz in steps of .2 Hz
            cfg.t_ftimwin    =  ones(length(cfg.foi),1).*10;   % length of time window =6 sec
            cfg.keeptrials   = 'yes';
            TFdataBefore     = ft_freqanalysis(cfg, dataBefore);
            TFdataAFTER      = ft_freqanalysis(cfg, dataAfter);
    
            nFc=nFc+1;
            av_PowDataEOBefore(nFc,:,:) = squeeze(10*log10(nanmean(TFdataBefore.powspctrm,4)));
            av_PowDataEOAfter(nFc,:,:) = squeeze(10*log10(nanmean(TFdataAfter.powspctrm,4)));
%     if redo==1 || exist([data_path filesep 'Preproc' filesep 'CIcf_ft_' file_name(1:end-4) '.mat'])==0
%                save([data_path filesep 'Preproc' filesep 'CIcf_ft_' file_name(1:end-4)],'data');
% 
%     end
    
end
    
plot(TFdataBefore.freq,squeeze(mean(av_PowDataEOBefore(:,match_str(data.label,'Oz'),:),1)))
plot(TFdataAfter.freq,squeeze(mean(av_PowDataEOAfter(:,match_str(data.label,'Oz'),:),1)))

