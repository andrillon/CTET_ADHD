%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

% path_SPM12='/Users/tand0009/Work/local/spm12/';
path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
path_localsleep='/Users/tand0009/WorkGit/projects/inprogress/wanderIM/localsleep';
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';

save_path='/Users/tand0009/Data/ADHD_CTET/';
root_path='/Volumes/ANDRILLON_HD1/ADHD_CTET/';
% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(path_SPM12))
addpath(genpath(path_LSCPtools))
addpath(path_localsleep)
addpath(path_fieldtrip)
ft_defaults;

% select relevant files, here baseline blocks
adhd_files=dir([root_path filesep 'adhds' filesep '*' filesep '*.bdf']);
controls_files=dir([root_path filesep 'controls' filesep '*' filesep '*.bdf']);

%% loop across trials for baseline blocks
redo=0;
for nF=1:length(adhd_files)
    % load file with spm
    if strcmp(adhd_files(nF).name,'CB_CTET.bdf')
        continue;
    end
    fprintf('... file: %s\n',adhd_files(nF).name)
    
    SubID=adhd_files(nF).name;
    SubID=SubID(1:findstr(SubID,'.')-1);
    if redo==0 && exist([save_path filesep 'SW_detection' filesep 'ADHD_CTET_SW_' SubID '.mat'])~=0
        continue;
    end
    hdr = ft_read_header([adhd_files(nF).folder filesep adhd_files(nF).name]);
    dat = ft_read_data([adhd_files(nF).folder filesep adhd_files(nF).name],'header',hdr);
    evt = ft_read_event([adhd_files(nF).folder filesep adhd_files(nF).name]);
    fprintf('... ... duration according to EEG data %g hours\n',size(dat,2)/hdr.Fs/60/60)
    
    numEpochs=floor(size(dat,2)/30/hdr.Fs);
    all_Waves=[];
    
    temp_data=dat(1:64,:);
    temp_data=temp_data-repmat(mean(temp_data(match_str(hdr.label,{'TP7','TP8'}),:),1),size(temp_data,1),1);
    temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
    
    [twa_results]=twalldetectnew_TA_v2(temp_data,hdr.Fs,0);
    for nE=1:64
        all_Waves=[all_Waves ; [repmat([1 nF nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
            cell2mat(twa_results.channels(nE).negzx)' ...
            cell2mat(twa_results.channels(nE).poszx)' ...
            cell2mat(twa_results.channels(nE).wvend)' ...
            cell2mat(twa_results.channels(nE).maxnegpk)' ...
            cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
            cell2mat(twa_results.channels(nE).maxpospk)' ...
            cell2mat(twa_results.channels(nE).maxpospkamp)' ...
            cell2mat(twa_results.channels(nE).mxdnslp)' ...
            cell2mat(twa_results.channels(nE).mxupslp)' ...
            cell2mat(twa_results.channels(nE).maxampwn)' ...
            cell2mat(twa_results.channels(nE).minampwn)' ...
            ]];
    end
    fprintf('\n')
    save([save_path filesep 'SW_detection' filesep 'ADHD_CTET_allSW_' SubID],'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=75/2;
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_detection' filesep 'ADHD_CTET_SW_' SubID],'slow_Waves','hdr','paramSW')
    
end

%% loop across trials for baseline blocks
for nF=1:length(controls_files)
    % load file with spm
    
    fprintf('... file: %s\n',controls_files(nF).name)
    
    SubID=controls_files(nF).name;
    SubID=SubID(1:findstr(SubID,'.')-1);
    if strcmp(controls_files(nF).name,'AK_EO_Before.bdf')
        continue;
    end
    if redo==0 && exist([save_path filesep 'SW_detection' filesep 'ADHD_CTET_SW_' SubID '.mat'])~=0
        continue;
    end
       hdr = ft_read_header([controls_files(nF).folder filesep controls_files(nF).name]);
 dat = ft_read_data([controls_files(nF).folder filesep controls_files(nF).name],'header',hdr);
    evt = ft_read_event([controls_files(nF).folder filesep controls_files(nF).name]);
    fprintf('... ... duration according to EEG data %g hours\n',size(dat,2)/hdr.Fs/60/60)
    
    numEpochs=floor(size(dat,2)/30/hdr.Fs);
    all_Waves=[];
    
    temp_data=dat(1:64,:);
    temp_data=temp_data-repmat(mean(temp_data(match_str(hdr.label,{'TP7','TP8'}),:),1),size(temp_data,1),1);
    temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
    
    [twa_results]=twalldetectnew_TA_v2(temp_data,hdr.Fs,0);
    for nE=1:64
        all_Waves=[all_Waves ; [repmat([2 nF nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
            cell2mat(twa_results.channels(nE).negzx)' ...
            cell2mat(twa_results.channels(nE).poszx)' ...
            cell2mat(twa_results.channels(nE).wvend)' ...
            cell2mat(twa_results.channels(nE).maxnegpk)' ...
            cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
            cell2mat(twa_results.channels(nE).maxpospk)' ...
            cell2mat(twa_results.channels(nE).maxpospkamp)' ...
            cell2mat(twa_results.channels(nE).mxdnslp)' ...
            cell2mat(twa_results.channels(nE).mxupslp)' ...
            cell2mat(twa_results.channels(nE).maxampwn)' ...
            cell2mat(twa_results.channels(nE).minampwn)' ...
            ]];
    end
    fprintf('\n')
    save([save_path filesep 'SW_detection' filesep 'ADHD_CTET_allSW_' SubID],'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=75/2;
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_detection' filesep 'ADHD_CTET_SW_' SubID],'slow_Waves','hdr','paramSW')
    
end
