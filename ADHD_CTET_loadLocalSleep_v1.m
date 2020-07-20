%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';

save_path='/Users/tand0009/Data/ADHD_CTET/';
root_path='/Volumes/ANDRILLON_HD1/ADHD_CTET/';
path_export='/Users/tand0009/Work/local/export_fig/';

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(path_LSCPtools))
addpath(genpath(path_export))

% select relevant files, here baseline blocks
list_files=dir([save_path filesep 'SW_detection' filesep 'ADHD_CTET_SW_*EO*.mat']);

path_fig='/Users/tand0009/Work/Documents/Grants/NHMRC/2020/Ideas_LocalSleep2/material';

%% loop across trials for baseline blocks
all_SW=[];
all_Counts=[];
thrAbs=50;
for nF=1:length(list_files)
%    if nF==11
%        continue;
%    end
    load([list_files(nF).folder filesep list_files(nF).name]); %,'slow_Waves','hdr','paramSW')
    slow_Waves(slow_Waves(:,4)<50,:)=[];
    all_SW=[all_SW ; slow_Waves];
    
    [nout,xout]=histcounts(slow_Waves(:,3),64);
    all_Counts=[all_Counts ; [slow_Waves(1,1:3) nout/(hdr.nSamples/hdr.Fs/60)]];
end


%%
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

cfg = [];
cfg.layout = 'biosemi64.lay';
layout=ft_prepare_layout(cfg);


limMax=5;


figure; 
temp_topo=nanmean(all_Counts(all_Counts(:,1)==1,4:65));
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title({'ADHD-OFF',sprintf('N=%g',sum(all_Counts(:,1)==1))}); h=colorbar;  ylabel(h, 'waves/min')
caxis([0 1]*limMax)
set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;
set(h,'FontSize',22);

% export_fig([path_fig filesep 'CTETADHD_ADHD_OFF_SWdensity.fig'])
% export_fig([path_fig filesep 'CTETADHD_ADHD_OFF_SWdensity.eps'],'-r 300')


figure; 
temp_topo=nanmean(all_Counts(all_Counts(:,1)==2,4:65));
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title({'CONTROLS',sprintf('N=%g',sum(all_Counts(:,1)==2))}); h=colorbar;  ylabel(h, 'waves/min')
caxis([0 1]*limMax)
set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;
set(h,'FontSize',22);

% export_fig([path_fig filesep 'CTETADHD_CONTROLS_SWdensity.fig'])
% export_fig([path_fig filesep 'CTETADHD_CONTROLS_SWdensity.eps'],'-r 300')


ADHD_CTET_allE{1}=nanmean(all_Counts(all_Counts(:,1)==1,4:65),2);
ADHD_CTET_allE{2}=nanmean(all_Counts(all_Counts(:,1)==2,4:65),2);
% save([path_fig filesep 'CTETADHD_res.mat'])


figure; 
[h,pV,~,stats]=ttest2(all_Counts(all_Counts(:,1)==1,4:65),all_Counts(all_Counts(:,1)==2,4:65));
temp_topo=stats.tstat;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title({'ADHD vs Controls','unpaired T-test'}); h=colorbar;  ylabel(h, 't-value')
caxis([-1 1]*4)
set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;
set(h,'FontSize',22);

% export_fig([path_fig filesep 'CTETADHD_CONTROLS_SWdensity.fig'])
% export_fig([path_fig filesep 'CTETADHD_CONTROLS_SWdensity.eps'],'-r 300')

