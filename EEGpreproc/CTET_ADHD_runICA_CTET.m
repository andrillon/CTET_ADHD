%%
clear all;
close all;

%%
run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;

files=dir([data_path filesep '*' filesep '*' filesep '*CTET*.bdf']);

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

badChannels_table=readtable('../ADHD_CTET_task_BadChannels.csv');

%% loop on subjects
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    % %     % skip if not before file
    % %     findbefore=findstr(file_name,'efore');
    % %     if isempty(findbefore)
    % %         continue;
    % %     end
    
    if exist([data_path filesep 'Preproc' filesep 'Icfe_ft_' file_name(1:end-4) '.mat'])==0
        
        hdr=ft_read_header([folder_name filesep file_name]);
        load([data_path filesep 'Preproc' filesep 'fe_ft_' file_name(1:end-4)]);
        
        %%% take out trial
        thisF=match_str(badChannels_table.FileName,['fe_ft_' file_name(1:end-4)]);
        tempChannels=badChannels_table.Bad_Channels{thisF}; tempChannels(tempChannels==' ')=[];
        eval(sprintf('badChannels={%s};',tempChannels));
        badChannels(cellfun(@isempty,badChannels))=[];
        eval(['badTrials=[' badChannels_table.Bad_Trials{thisF} '];']);
        
        cfg=[];
        cfg.trials          = setdiff(1:length(data.trial),badTrials);
        data = ft_preprocessing(cfg, data);
        
        if ~isempty(badChannels)
            fprintf('... ... interpolating %g channels\n',length(badChannels))
            % find neighbours
            cfg=[];
            cfg.method        = 'triangulation';
            cfg.layout        = layout;
            cfg.feedback      = 'no';
            cfg.channel = layout.label;
            [neighbours] = ft_prepare_neighbours(cfg);
            
            % interpolate channels
            cfg=[];
            cfg.method         = 'weighted';
            cfg.badchannel     = badChannels;
            cfg.missingchannel = [];
            cfg.neighbours     = neighbours;
            cfg.trials         = 'all';
            cfg.layout         = layout;
            cfg.channel = layout.label;
            [data] = ft_channelrepair(cfg, data);
        end
        
        cfg=[];
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        data = ft_preprocessing(cfg,data);
        
        %%% run ICA
        rankICA = rank(data.trial{1,1});
        cfg        = [];
        cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
        cfg.numcomponent = rankICA;
        comp = ft_componentanalysis(cfg, data);
        save([data_path filesep 'Preproc' filesep 'Icfe_ft_' file_name(1:end-4)],'data','comp','rankICA');
    end
    toc;
end