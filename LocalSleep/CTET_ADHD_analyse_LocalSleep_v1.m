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
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

ChanLabels=layout.label(1:64);
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
    
%     load([data_path filesep 'Preproc' filesep 'CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    load([data_path filesep 'Preproc' filesep 'fixThr_CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
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
        sub_table_behav=table_behav(table_behav.BlockN==nBl,:);
        duration_block=(sub_table_behav.Sample(end)-sub_table_behav.Sample(1))/hdr.Fs/60;
        densSW=nout/duration_block;
        all_slowWaves(nFc,nBl,:)=densSW;
        
        these_Amplitude=slow_Waves(slow_Waves(:,2)==nBl,4);
        byElec_Amplitude=nan(1,64);
        for nE=1:64
            byElec_Amplitude(nE)=mean(these_Amplitude(these_SWelectrodes==nE));
        end
         all_slowWaves_P2P(nFc,nBl,:)=byElec_Amplitude;
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
% 1: Split between ADHD and Controls

% 2: For Electrode Cz (48), plot the average SW density across blocks
% (simpleTplot or RainCloudPlot) squeeze(all_slowWaves(:,:,match_str(ChanLabels,'Cz')))

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves(:,:,match_str(ChanLabels,'Cz'))),0,'k',0,'-',0.5,1,0,1,2);
title('Averaged slow-waves density across blocks at Cz');

% 2b: Do the same thing for ADHD and controls separately

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Cz'))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Cz'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('Averaged slow-waves density across blocks at Cz');

% 3: Plot scalp topographies of average SW density (simpleTopoPlot_ft)
% squeeze(mean(mean(all_slowWaves,1),2))

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

SW_topo=squeeze(nanmean(nanmean(all_slowWaves,1),2));
figure;
subplot(1,3,1);
simpleTopoPlot_ft(SW_topo', layout,'labels',[],0,1);
colorbar;
title('Topography of average SW density')
caxis([0 4.5])
% 3b: Do the same for ADHD and controls separately

SW_topo=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'Control'),:,:))));
subplot(1,3,2);
simpleTopoPlot_ft(SW_topo', layout,'labels',[],0,1);
colorbar
title('Topography of average SW density for Controls')
caxis([0 4.5])

SW_topo=squeeze(nanmean(nanmean(all_slowWaves(match_str(group_SW,'ADHD'),:,:))));
subplot(1,3,3);
simpleTopoPlot_ft(SW_topo', layout,'labels',[],0,1);
colorbar;
title('Topography of average SW density for ADHDs')
caxis([0 4.5])

%%
figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(:,:,match_str(ChanLabels,'Cz'))),0,'k',0,'-',0.5,1,0,1,2);
title('Averaged slow-waves P2P across blocks at Cz');

% 2b: Do the same thing for ADHD and controls separately

figure;
hp=[];
[~,hp(1)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(match_str(group_SW,'Control'),:,match_str(ChanLabels,'Cz'))),0,'b',0,'-',0.5,1,0,1,2);
hold on;
[~,hp(2)]=simpleTplot(1:8,squeeze(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,match_str(ChanLabels,'Cz'))),0,'r',0,'-',0.5,1,0,1,2);
hold on;
legend(hp,{'Controls','ADHDs'})
title('Averaged slow-waves P2P across blocks at Cz');

%%
SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P,1),2));
figure;
subplot(1,3,1);
simpleTopoPlot_ft(SW_topo, layout,'labels',[],0,1);
colorbar;
title('Topography of average SW P2P')
% caxis([20 42])
% 3b: Do the same for ADHD and controls separately

SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'Control'),:,:))));
subplot(1,3,2);
simpleTopoPlot_ft(SW_topo, layout,'labels',[],0,1);
colorbar
title('Topography of average SW P2P for Controls')
% caxis([20 42])

SW_topo=squeeze(nanmean(nanmean(all_slowWaves_P2P(match_str(group_SW,'ADHD'),:,:))));
subplot(1,3,3);
simpleTopoPlot_ft(SW_topo, layout,'labels',[],0,1);
colorbar;
title('Topography of average SW P2P for ADHDs')
% caxis([20 42])
