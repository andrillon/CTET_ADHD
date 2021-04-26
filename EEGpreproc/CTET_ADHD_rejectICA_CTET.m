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

load('../ADHD_CTET_Task_ICA.mat');

%% loop on subjects
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
    
    if exist([data_path filesep 'Preproc' filesep 'CIcfe_ft_' file_name(1:end-4) '.mat'])==0
        
        hdr=ft_read_header([folder_name filesep file_name]);
        load([data_path filesep 'Preproc' filesep 'Icfe_ft_' file_name(1:end-4)]);
        fprintf('... processing %s\n',file_name(1:end-4));
        data_ori=data;
        
        %%% Reject bad component
        cfg = [];
        this_line=match_str(ICA_table.FileName,['Icfe_ft_' file_name(1:end-4)]);
        if isempty(this_line)
            warning(sprintf('... missing %s in ICA detection\n',file_name(1:end-4)));
            continue;
        end
        these_components=ICA_table.blink{this_line};
        these_components2=ICA_table.saccade{this_line};
        if isempty(these_components2)
            eval(sprintf('cfg.component = [%s];',these_components)); % to be removed component(s)
        else
            eval(sprintf('cfg.component = [%s,%s];',these_components,these_components2)); % to be removed component(s)
        end
        data = ft_rejectcomponent(cfg, comp, data);
        
        cfg=[];
        cfg.reref           = 'yes';
        cfg.refchannel      = 'all';
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.2 0];
        data = ft_preprocessing(cfg,data);
        
        save([data_path filesep 'Preproc' filesep 'CIcfe_ft_' file_name(1:end-4)],'data');
    end
    toc;
end