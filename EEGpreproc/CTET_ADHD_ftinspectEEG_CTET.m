%%
run ../localdef_ADHD_CTET.m
addpath((path_fieldtrip));
ft_defaults;

%% choose and load subject
List_Subj=dir([data_path filesep 'Preproc' filesep 'fe_ft_*.mat']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
fprintf('... working on %s\n',ListNames{pick})
load([data_path filesep 'Preproc' filesep ListNames{pick}])
 
%% reject trials
cfg          = [];
cfg.method   = 'summary';
cfg.alim     = 5e-5;
data2        = ft_rejectvisual(cfg,data);

%% display data
cfg=[];
cfg.continuous='no';
cfg.allowoverlap='true';
cfg.viewmode='vertical';
cfg = ft_databrowser(cfg, data);
