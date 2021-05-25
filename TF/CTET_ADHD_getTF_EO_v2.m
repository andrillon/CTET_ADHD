%%
clear all;
close all;

%%
run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools));
files=dir([data_path filesep 'Preproc' filesep 'CIcf_ft_*.mat']);

f_range = [2, 30];
settings = struct();  % Use defaults
addpath(genpath(fooof_path));

%% loop on subjects
redo=0;
nFc=0;
all_SubIDs=[];
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubIDlong=file_name(9:end-4);
    SubID=file_name(9:end-4);
    seps=findstr(SubID,'_');
    if length(seps)>1
    SubID=SubID(1:seps(end-1)-1);
    else
    SubID=SubID(1:seps(1)-1);
    end
    tic;
    
    if isempty(findstr(file_name,'efore'))==0 || strcmp(file_name,'CIcf_ft_FA_EO.mat')
    nFc=nFc+1;
        %         av_PowDataEO_Before(nFc1,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
        all_SubIDs{nFc}=SubID;
    else
        continue;
    end
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
    
    av_PowDataEO(nFc,1,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
    
    temp_Power=squeeze((nanmean(TFdata.powspctrm,4)));
    %temp_fooof=[];
    %for nE=1:64
        %fooof_results = fooof(TFdata.freq, temp_Power(nE,:), f_range, settings,1);
        %temp_fooof(nE,:)=fooof_results.background_params;
    %end
    %av_PowDataEO_FOOF(nFc,:,:) =temp_fooof;
    
    
    
    orifoldername=orifile.folder;
    if isempty(findstr(orifoldername,'controls'))==0
        group_PowDataEO{nFc}='Control';
        design_PowDataEO(1,nFc)=0;
    elseif isempty(findstr(orifoldername,'adhds'))==0
        group_PowDataEO{nFc}='ADHD';
        design_PowDataEO(1,nFc)=1;
    end
end
%%
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubIDlong=file_name(9:end-4);
    SubID=file_name(9:end-4);
    seps=findstr(SubID,'_');
if length(seps)>1
    SubID=SubID(1:seps(end-1)-1);
    else
    SubID=SubID(1:seps(1)-1);
end
tic;
    
         thisF=find(ismember(all_SubIDs,SubID));
   if isempty(findstr(file_name,'fter'))==0 && isempty(thisF)==0
        %         av_PowDataEO_Before(nFc1,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
    else
        continue;
    end
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
    
    av_PowDataEO(thisF,2,:,:) = squeeze(10*log10(nanmean(TFdata.powspctrm,4)));
 
end

%% Matching electrodes
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
        matching_elec(nE)=(match_str(TFdata.label,layout.label(nE)));
    end

%Plot topo report
temp_topo=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),2,matching_elec,TFdata.freq>4 & TFdata.freq<7),4),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colorbar;
title('Topography on the [4 7]Hz frequency domain for Controls after the CTET task')
caxis([-13 -5])

%Plot t-values topo
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:size(av_PowDataEO,3)
    A=squeeze(nanmean(av_PowDataEO(match_str(group_PowDataEO,'Control'),1,nE,TFdata.freq>4 & TFdata.freq<7),4));
    B=squeeze(nanmean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),1,nE,TFdata.freq>4 & TFdata.freq<7),4));
   [h,pV,~,stat]=ttest2(B,A); 
   temp_topo_tval(nE)=stat.tstat;
   temp_topo_pval(nE)=pV;
end
figure;
simpleTopoPlot_ft(temp_topo_tval(matching_elec),layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Topography difference TF Power before CTET ADHD/Control (tvalue)')
caxis([-1 1]*4)

%% 1: Plot average power for Cz Oz and Fz
% figure;
% plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Fz'),:),1)))
% hold on;
% plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Cz'),:),1)))
% hold on;
% plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Pz'),:),1)))
% hold on;
% plot(TFdata.freq,squeeze(mean(av_PowDataEO(:,match_str(TFdata.label,'Oz'),:),1)))
% hold on;
% legend({'Fz','Cz','Pz', 'Oz'});
% title('Average power spectrum : eyes open');
% xlabel('Frequency (Hz)')
% ylabel('dB Power')

%% 2: Split between before and after
figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(:,1,match_str(TFdata.label,'Fz'),:))),0,[1 1 1]*0.5,0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(:,2,match_str(TFdata.label,'Fz'),:))),0,[1 1 1]*0,0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Before CTET','After CTET'})
title('Average power spectrum at electrode Oz: eyes open, before vs after the CTET task');
xlabel('Frequency (Hz)')
ylabel('dB Power')
format_fig;

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(:,2,match_str(TFdata.label,'Fz'),:)))-...
squeeze((av_PowDataEO(:,1,match_str(TFdata.label,'Fz'),:))),0,[1 1 1]*0.5,[2 0.05 0.05 1000],'-',0.5,1,0,1,2);

%% 3: Split between ADHD and controls

Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'Control'),2,match_str(TFdata.label,'Fz'),:))),0,Colors(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'ADHD'),2,match_str(TFdata.label,'Fz'),:))),0,Colors(2,:),0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHD'})
title('Average power spectrum at electrode Fz after the CTET task');
xlabel('Frequency (Hz)')
ylabel('dB Power')
%% 4: Split between before and after + ADHD and controls

figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'Control'),1,match_str(TFdata.label,'Fz'),:))),0,Colors(1,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'ADHD'),1,match_str(TFdata.label,'Fz'),:))),0,Colors(2,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(3)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'Control'),2,match_str(TFdata.label,'Fz'),:))),0,Colors(3,:),0,'-',0.5,1,0,1,2);
hold on;
[~,hp(4)]=simpleTplot(TFdata.freq,squeeze((av_PowDataEO(match_str(group_PowDataEO,'ADHD'),2,match_str(TFdata.label,'Fz'),:))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls before','ADHD before','Controls after', 'ADHD after'})
title('Average power spectrum at electrode Fz after the CTET task');
xlabel('Frequency (Hz)')
ylabel('dB Power')


%%
temp_topo=squeeze(mean(mean(av_PowDataEO(:,matching_elec,TFdata.freq>8.5 & TFdata.freq<11.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colorbar;
title('Alpha [8.5 11.5]')


temp_topo=squeeze(mean(mean(av_PowDataEO(:,matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colorbar;
title('Theta [4.5 7.5]')


temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo=temp_topo_C-temp_topo_A;
figure;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast Controls/ADHDs')

%Theta contrast ADHDs > Controls
temp_topo_C=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'Control'),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo_A=squeeze(mean(mean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo2=temp_topo_A-temp_topo_C;
figure;
simpleTopoPlot_ft(temp_topo2', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast ADHDs/Controls')

%Comparisons

% cfg = [];
% cfg.layout = 'biosemi64.lay';
% cfg.channel  = TFdata.label;
% cfg.center      = 'yes';
% layout=ft_prepare_layout(cfg);

%Theta Controls Before
temp_topoCB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'Before')),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topoCB', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] Controls Before')

%Theta Controls After
temp_topoCA=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'After')),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
figure;
simpleTopoPlot_ft(temp_topoCA', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] Controls After')

%Theta contrast Controls After > Before
temp_topoCB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'Before')),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topoCA=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'Controls'),match_str(cond_PowDataEO,'After')),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
temp_topo3=temp_topoCA-temp_topoCB;
figure;
simpleTopoPlot_ft(temp_topo3', layout,'labels',[],0,1);
colorbar;
title('Theta [4.5 7.5] contrast Controls After/Before')


%Theta ADHDs Before
temp_topoAB=squeeze(mean(mean(av_PowDataEO(intersect(match_str(group_PowDataEO,'ADHD'),match_str(cond_PowDataEO,'Before')),matching_elec,TFdata.freq>4.5 & TFdata.freq<7.5),3),1));
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
title('Theta [4.5 7.5]')

%% Permutation
% cfg = [];
% cfg.channel          = 'all';
% cfg.latency          = 'all';
% cfg.frequency        = 'all';
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_indepsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 500;
% % prepare_neighbours determines what sensors may form clusters
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, dataFC);
% 
% design = zeros(1,size(freqFIC_planar_cmb.powspctrm,1) + size(freqFC_planar_cmb.powspctrm,1));
% design(1,1:size(freqFIC_planar_cmb.powspctrm,1)) = 1;
% design(1,(size(freqFIC_planar_cmb.powspctrm,1)+1):(size(freqFIC_planar_cmb.powspctrm,1)+...
% size(freqFC_planar_cmb.powspctrm,1))) = 2;
% 
% cfg.design           = design;
% cfg.ivar             = 1;
% 
% [stat] = ft_freqstatistics(cfg, freqFIC_planar_cmb, freqFC_planar_cmb);

%% Fig for report
% Difference Before/After
Bef = squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),match_str(TFdata.label,'Fz'),:),1))
Af = squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),match_str(TFdata.label,'Fz'),:),1))
diffAB = Af - Bef
figure;
plot(TFdata.freq,diffAB)

%retester ça avec nanmean
Bef2 = squeeze(nanmean(av_PowDataEO(match_str(group_PowDataEO,'Control'),match_str(TFdata.label,'Fz'),:)))
Af2 = squeeze(nanmean(av_PowDataEO(match_str(group_PowDataEO,'ADHD'),match_str(TFdata.label,'Fz'),:)))
diff2AB = Af2 - Bef2
figure;
hp=[];
[~,hp(1)]=simpleTplot(TFdata.freq,diff2AB,0,'r',2,'-',0.5,1,0,1,2);

%[h, pV_diffGroup,~,stat_diffGroup]=ttest2(Af,Bef);
%fprintf('... unpaired t-test between cond Before/After on Fz: p=%g, t-value=%g, df=%g\n',...
    %pV_diffGroup,stat_diffGroup.tstat,stat_diffGroup.df)
    
% cfg = [];
% cfg.layout = 'biosemi64.lay';
% cfg.channel  = TFdata.label;
% cfg.center      = 'yes';
% layout=ft_prepare_layout(cfg);

%cfg = [];
%cfg.layout = 'biosemi64.lay';
%cfg.channel = layout.label;
%cfg.channel(match_str(layout.label,{'Iz','P7','P8'}))=[];
%cfg.center      = 'yes';
%layout=ft_prepare_layout(cfg);

temp_topo_Bef=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'Before'),:,TFdata.freq>8 & TFdata.freq<10),3),1));
temp_topo_Af=squeeze(mean(mean(av_PowDataEO(match_str(cond_PowDataEO,'After'),:,TFdata.freq>8 & TFdata.freq<10),3),1));
temp_topo2=temp_topo_Af-temp_topo_Bef;

figure;
simpleTopoPlot_ft(temp_topo2', layout,'on',[],0,1);
colorbar;
title('[8 10] contrast After/Before')

figure;
simpleTopoPlot_ft(temp_topo_Bef', layout,'on',[],0,1);
colorbar;
title('[8 10] Before')
caxis([-11 -3])

figure;
simpleTopoPlot_ft(temp_topo_Af', layout,'on',[],0,1);
colorbar;
title('[8 10] After')
caxis([-11 -3])

%%
% Fig for report Control/ADHD
new = av_PowDataEO(match_str(group_PowDataEO,'ADHD'),:,TFdata.freq>8 & TFdata.freq<10)-squeeze(mean(av_PowDataEO(match_str(cond_PowDataEO,'Control'),:,TFdata.freq>8 & TFdata.freq<10)))