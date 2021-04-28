%%
clear all;
close all;

%%
run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;

files=dir([data_path filesep 'Preproc' filesep 'CIcf_ft_*.mat']);

f_range = [2, 30];
settings = struct();  % Use defaults
addpath(genpath(fooof_path));

%% loop on subjects
redo=0;
nFc=0;
nFc1=0;
nFc2=0;
nFc3=0;
nFc4=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubIDlong=file_name(9:end-4);
    SubID=file_name(9:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    if redo==1 || exist([data_path filesep 'Preproc' filesep 'TF_' file_name])==0
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
        cfg.toi         =  1:5:data.time{1}(end);                         % analysis 2 to 30 Hz in steps of .2 Hz
        cfg.t_ftimwin    =  ones(length(cfg.foi),1).*10;   % length of time window =6 sec
        cfg.keeptrials   = 'yes';
        TFdata           = ft_freqanalysis(cfg, data);
        
        save([data_path filesep 'Preproc' filesep 'TF_' file_name(1:end-4)],'TFdata');
    else
             orifile=dir([data_path filesep '*' filesep  '*' filesep SubIDlong '.bdf']);
   load([data_path filesep 'Preproc' filesep 'TF_' file_name(1:end-4)]);
    end
    
    nFc=nFc+1;
    av_PowDataEO(nFc,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
    
    temp_Power=squeeze((nanmean(TFdata.powspctrm,4)));
    temp_fooof=[];
    for nE=1:64
        fooof_results = fooof(TFdata.freq, temp_Power(nE,:), f_range, settings,1);
        temp_fooof(nE,:)=fooof_results.background_params;
    end
    av_PowDataEO_FOOF(nFc,:,:) =temp_fooof;
    
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
    if isempty(findstr(orifoldername,'controls'))==0
        nFc3=nFc3+1;
        group_PowDataEO{nFc}='Control';
    elseif isempty(findstr(orifoldername,'adhds'))==0
        nFc4=nFc4+1;
        group_PowDataEO{nFc}='ADHD';
    end
end

%% 1: Plot average power for Cz Oz and Fz
figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Fz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Cz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Pz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(data.label,'Oz'),:),1)))
hold on;
legend({'Fz','Cz','Pz', 'Oz'});
title('Average power spectrum : eyes open');
xlabel('Frequency (Hz)')
ylabel('dB Power')

%% 2: Split between before and after
figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),match_str(data.label,'Oz'),:),1)),'Color','b')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),match_str(data.label,'Oz'),:),1)),'Color','r')
hold on;
legend({'Before CTET','After CTET'})
title('Average power spectrum at electrode Oz: eyes open, before vs after the CTET task');
xlabel('Frequency (Hz)')
ylabel('dB Power')

%% 3: Split between ADHD and controls

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'Control'),match_str(data.label,'Fz'),:))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'ADHD'),match_str(data.label,'Fz'),:))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHD'})
title('Average power spectrum at electrode Fz: eyes open, controls vs ADHDs');
xlabel('Frequency (Hz)')
ylabel('dB Power')
%% 4: Split between before and after + ADHD and controls

figure;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Control'),match_str(cond_PowDataEO,'Before')),match_str(data.label,'Fz'),:),1)),'Color','g')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Control'),match_str(cond_PowDataEO,'After')),match_str(data.label,'Fz'),:),1)),'Color','c')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'Before')),match_str(data.label,'Fz'),:),1)),'Color','b')
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'After')),match_str(data.label,'Fz'),:),1)),'Color','m')
hold on;
legend({'Controls Before','Controls After','ADHD Before','ADHD After'})
title('Average power spectrum at electrode Fz: eyes open');
xlabel('Frequency (Hz)')
ylabel('dB Power')

%%

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel  = TFdata.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);


temp_topo=squeeze(mean(mean(av_PowDataEO(:,:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5]')


temp_topo=squeeze(mean(mean(av_PowDataEO(:,:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5]')


temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo=temp_topo_C-temp_topo_A;
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5]')