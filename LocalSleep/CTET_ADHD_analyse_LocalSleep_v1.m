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
    
    load([data_path filesep 'Preproc' filesep 'CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    % 1: Subject Number
    % 2: Block Number
    % 3: Electrode Number
    % 4: P2P amplitude
    % 5: Start slow wave (sample from block onset)
    nFc=nFc+1;
    for nBl=1:8
        these_SWelectrodes=slow_Waves(slow_Waves(:,2)==nBl,3);
        nout=hist(these_SWelectrodes,1:64);
        sub_table_behav=table_behav(table_behav.BlockN==nBl,:);
        duration_block=(sub_table_behav.Sample(end)-sub_table_behav.Sample(1))/hdr.Fs/60;
        densSW=nout/duration_block;
        all_slowWaves(nFc,nBl,:)=densSW;
    end
end

%%
% 1: Split between ADHD and Controls

% 2: For Electrode Cz (48), plot the average SW density across blocks
% (simpleTplot or RainCloudPlot) squeeze(all_slowWaves(:,:,match_str(ChanLabels,'Cz')))

% 2b: Do the same thing for ADHD and controls separately

% 3: Plot scalp topographies of average SW density (simpleTopoPlot_ft)
% squeeze(mean(mean(all_slowWaves,1),2))

% 3b: Do the same for ADHD and controls separately