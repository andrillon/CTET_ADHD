%%
clear all;
close all;

%%
run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;

files=dir([data_path filesep '*' filesep '*' filesep '*EO*.bdf']);

%% loop on subjects
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'_');
    SubID=SubID(1:seps(1)-1);
    tic;
    fprintf('... working on %s (%g/%g)\n',file_name,nF,length(files))
    
    if exist([data_path filesep 'Preproc' filesep 'f_ft_' file_name(1:end-4) '.mat'])==0
        
        hdr=ft_read_header([folder_name filesep file_name]);
        
        %%% Define epochs
        cfg=[];
        cfg.SubID               = SubID;
        cfg.dataset             = [folder_name filesep file_name];
        
        cfg.channel        = hdr.label(match_str(hdr.chantype,'eeg'));
        cfg.demean         = 'yes';
        cfg.lpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord      = 4;
        cfg.lpfreq         = 40;
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.hpfilttype     = 'but';
        cfg.hpfiltord      = 4;
        cfg.hpfreq         = 0.1;
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50 100 150]; % set up the frequencies for notch filtering
        
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        
        dat                = ft_preprocessing(cfg); % read raw data
        
        cfgbs=[];
        cfgbs.resamplefs      = 256;
        cfgbs.detrend         = 'no';
        cfgbs.demean          = 'yes';
        cfgbs.baselinewindow  = [-0.2 0];
        data                  = ft_resampledata(cfgbs,dat); % read raw data
        save([data_path filesep 'Preproc' filesep 'f_ft_' file_name(1:end-4)],'data');
        
        %         cfg2=[];
        %         av_data = ft_timelockanalysis(cfg2, data);
        %         cfg3=[];
        %         cfg3.channel = 'Cz';
        %         figure; ft_singleplotER(cfg3,av_data);
    end
    toc;
end