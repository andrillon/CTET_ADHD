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
    %temp_fooof=[];
    %for nE=1:64
        %fooof_results = fooof(TFdata.freq, temp_Power(nE,:), f_range, settings,1);
        %temp_fooof(nE,:)=fooof_results.background_params;
    %end
    %av_PowDataEO_FOOF(nFc,:,:) =temp_fooof;
    
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
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Fz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Cz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Pz'),:),1)))
hold on;
plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Oz'),:),1)))
hold on;
legend({'Fz','Cz','Pz', 'Oz'});
title('Average power spectrum : eyes open');
xlabel('Frequency (Hz)')
ylabel('dB Power')

%% 2: Split between before and after

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(cond_PowDataEO,'Before'),match_str(TFdata.label,'Oz'),:))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(cond_PowDataEO,'After'),match_str(TFdata.label,'Oz'),:))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Before','After'})
title('Power spectrum at electrode Oz : eyes open, before vs after CTET task');
xlabel('Frequency (Hz)')
ylabel('dB Power')

%% 3: Split between ADHD and controls

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'Control'),match_str(TFdata.label,'Fp1'),:))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'ADHD'),match_str(TFdata.label,'Fp1'),:))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHD'})
title('Average power spectrum at electrode Fp1: eyes open, controls vs ADHDs');
xlabel('Frequency (Hz)')
ylabel('dB Power')
%% 4: Split between before and after + ADHD and controls

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(intersect(match_str(group_PowDataEO,'Control'),match_str(cond_PowDataEO,'Before')),match_str(TFdata.label,'Oz'),:))),0,'g',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(intersect(match_str(group_PowDataEO,'Control'),match_str(cond_PowDataEO,'After')),match_str(TFdata.label,'Oz'),:))),0,'c',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'Before')),match_str(TFdata.label,'Oz'),:))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'After')),match_str(TFdata.label,'Oz'),:))),0,'m',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls Before','Controls After','ADHD Before','ADHD After'})
title('Power spectrums at electrode Oz: eyes open');
xlabel('Frequency (Hz)')
ylabel('dB Power')

%% Topographies

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel  = TFdata.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

%Alpha topo
temp_topo=squeeze(mean(mean(av_PowDataEO(:,:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5]')

%Alpha before CTET
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] before CTET')

%Alpha after CTET
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] after CTET')

%Alpha contrast Before > After
temp_topo_B=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo_Af=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo=temp_topo_B-temp_topo_Af;
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] contrast Before/After')

%Alpha contrast After > Before
temp_topo_B=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo_Af=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo2=temp_topo_Af-temp_topo_B;
figure;
simpleTopoPlot_ft(temp_topo2', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5]contrast After/Before')

%Alpha Controls
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] Controls')

%Alpha ADHDS
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] ADHDS')

%Alpha contrast Controls > ADHDs
temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo=temp_topo_C-temp_topo_A;
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] contrast Controls/ADHDs')

%Alpha contrast ADHDs > Controls
temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
temp_topo2=temp_topo_A-temp_topo_C;
figure;
simpleTopoPlot_ft(temp_topo2', layout,'labels',[],0,1);
colorbar;
title('Alpha [8.5 11.5] contrast ADHDs/Controls')

%Theta topo
temp_topo=squeeze(mean(mean(av_PowDataEO(:,:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5]')

%Theta before CTET
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] before CTET')

%Theta after CTET
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] after CTET')

%Theta contrast Before > After
temp_topo_B=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_Af=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo=temp_topo_B-temp_topo_Af;
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast Before/After')

%Theta contrast After > Before
temp_topo_B=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_Af=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo2=temp_topo_Af-temp_topo_B;
figure;
simpleTopoPlot_ft(temp_topo2', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast After/Before')

%Theta Controls
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] Controls')

%Theta ADHDS
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] ADHDS')

%Theta contrast Controls > ADHDs
temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo=temp_topo_C-temp_topo_A;
figure;
simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast Controls/ADHDs')

%Theta contrast ADHDs > Controls
temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo2=temp_topo_A-temp_topo_C;
figure;
simpleTopoPlot_ft(temp_topo2', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast ADHDs/Controls')

%% Comparisons

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel  = TFdata.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

%Theta Controls Before
temp_topoCB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'Before')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topoCB', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] Controls Before')

%Theta Controls After
temp_topoCA=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'After')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topoCA', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] Controls After')

%Theta contrast Controls After > Before
temp_topoCB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'Before')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topoCA=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'After')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo3=temp_topoCA-temp_topoCB;
figure;
simpleTopoPlot_ft(temp_topo3', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast Controls After/Before')


%Theta ADHDs Before
temp_topoAB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'Before')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topoAB', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] ADHDs Before')

%Theta ADHDs After
temp_topoAA=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'After')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topoAA', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] ADHDs after')

%Theta contrast ADHDs After > Before
temp_topoAB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'Before')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topoAA=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'After')),:,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo4=temp_topoAA-temp_topoAB;
figure;
simpleTopoPlot_ft(temp_topo4', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast ADHDs After/Before')

%Theta contrast ADHDs After/Before // Controls After/Before

temp_topo5=temp_topo4-temp_topo3;
figure;
simpleTopoPlot_ft(temp_topo2', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast ADHD/Controls After/Before')