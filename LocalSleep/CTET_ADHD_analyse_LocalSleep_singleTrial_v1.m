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
layout=ft_prepare_layout(cfg);

cfg = []
cfg.layout = 'biosemi64.lay';
cfg.channel = layout.label;
cfg.channel(match_str(layout.label,{'Iz','P7','P8'}))=[];
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

ChanLabels=layout.label(1:end-2);
%%
all_table_behav_SW=[];
nFc=0;
group_SW=[];
table_SWandBehav=[];
table_P2PandBehav=[];
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
    load([data_path filesep 'Preproc' filesep  'CIcfeblock_ft_' file_name(1:end-4) '.mat']);
    matching_elec=[];
    for nE=1:length(layout.label)-2
        matching_elec(nE)=(match_str(data.label,layout.label(nE)));
    end
    
    %load([data_path filesep 'Preproc' filesep 'CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    load([data_path filesep 'Preproc' filesep 'relThrCTR_CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    %load([data_path filesep 'Preproc' filesep 'fixThr_CIcfeblock_ft_SW_' file_name(1:end-4)]); %,'slow_Waves','paramSW')
    % 1: Subject Number
    % 2: Block Number
    % 3: Electrode Number
    % 4: P2P amplitude
    % 5: Start slow wave (sample from block onset)
    % 12: Downward Slope
    % 13: Upward Slope
    nFc=nFc+1;
    table_behav_SW=[];
    for nBl=1:8
        sub_table_behav=table_behav(table_behav.BlockN==nBl,:);
        sub_table_behav.SampleBlock=sub_table_behav.Sample-sub_table_behav.Sample(1)+data.hdr.Fs;
        sub_table_behav.SampleBlockNext=[sub_table_behav.SampleBlock(2:end) ; NaN];
        for nE=1:64
            sub_table_behav.(data.label{nE})=zeros(size(sub_table_behav,1),1);
            sub_slow_Waves=slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,:);
            sub_slow_Waves(:,[5 6 7 8 10])=sub_slow_Waves(:,[5 6 7 8 10])*2;
            for nW=1:size(sub_slow_Waves,1)
                findTr=find(sub_slow_Waves(nW,5)>=sub_table_behav.SampleBlock & sub_slow_Waves(nW,5)<sub_table_behav.SampleBlockNext);
                if ~isempty(findTr)
                    sub_table_behav.(data.label{nE})(findTr)=1;
                end
            end
        end
        sub_table_behav.allE=sum(table2array(sub_table_behav(:,12:75)),2);
        table_behav_SW=[table_behav_SW ; sub_table_behav];
    end
    
    
    orifoldername=files(nF).folder;
    if isempty(findstr(orifoldername,'controls'))==0
        group_SW{nFc}='Control';
    elseif isempty(findstr(orifoldername,'adhds'))==0
        group_SW{nFc}='ADHD';
    end
    table_behav_SW.Group=cell(size(table_behav_SW,1),1);
    table_behav_SW.Group=repmat(group_SW(nFc),size(table_behav_SW,1),1);
    
    table_behav_SW.SubID=cell(size(table_behav_SW,1),1);
    table_behav_SW.SubID=repmat({SubID},size(table_behav_SW,1),1);
    
    writetable(table_behav_SW,[save_path filesep 'CTET_ADHD_behav_SW_' file_name(1:end-4) '.txt']);
    all_table_behav_SW=[all_table_behav_SW ; table_behav_SW];
end

%%
all_table_behav_SW.FA=1-all_table_behav_SW.corrTG;
all_table_behav_SW.Miss=1-all_table_behav_SW.corrNT;
clear coeff_*
fprintf('%2.0f/%2.0f\n',nE,length(data.label))
for nE=1:length(data.label)
    fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nE,length(data.label))
    mdl_FA=fitlme(all_table_behav_SW(~isnan(all_table_behav_SW.FA),:),sprintf('FA~1+BlockN+%s+(1+%s|SubID)',data.label{nE},data.label{nE}));
    coeff_FA(nE,1)=mdl_FA.Coefficients.tStat(match_str(mdl_FA.CoefficientNames,data.label{nE}));
    coeff_FA(nE,2)=mdl_FA.Coefficients.pValue(match_str(mdl_FA.CoefficientNames,data.label{nE}));
    
    mdl_Miss=fitlme(all_table_behav_SW(~isnan(all_table_behav_SW.Miss),:),sprintf('Miss~1+BlockN+%s+(1+%s|SubID)',data.label{nE},data.label{nE}));
    coeff_Miss(nE,1)=mdl_Miss.Coefficients.tStat(match_str(mdl_Miss.CoefficientNames,data.label{nE}));
    coeff_Miss(nE,2)=mdl_Miss.Coefficients.pValue(match_str(mdl_Miss.CoefficientNames,data.label{nE}));
    
    mdl_RT=fitlme(all_table_behav_SW(~isnan(all_table_behav_SW.RT),:),sprintf('RT~1+BlockN+%s+(1+%s|SubID)',data.label{nE},data.label{nE}));
    coeff_RT(nE,1)=mdl_RT.Coefficients.tStat(match_str(mdl_RT.CoefficientNames,data.label{nE}));
    coeff_RT(nE,2)=mdl_RT.Coefficients.pValue(match_str(mdl_RT.CoefficientNames,data.label{nE}));
end


%%
%T-values P2P amplitude topo
cmap_ttest=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap_ttest=flipud(cmap_ttest);

figure;
simpleTopoPlot_ft(coeff_FA(matching_elec,1),layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Effect of SW on FA')
caxis([-1 1]*5)
if ~isempty(find(coeff_FA(matching_elec,2)<fdr(coeff_FA(matching_elec,2),0.05)))
    ft_plot_lay_me(layout, 'chanindx',find(coeff_FA(matching_elec,2)<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

figure;
simpleTopoPlot_ft(coeff_Miss(matching_elec,1),layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Effect of SW on Miss')
caxis([-1 1]*5)
if ~isempty(find(coeff_Miss(matching_elec,2)<fdr(coeff_Miss(matching_elec,2),0.05)))
    ft_plot_lay_me(layout, 'chanindx',find(coeff_Miss(matching_elec,2)<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end

figure;
simpleTopoPlot_ft(coeff_RT(matching_elec,1),layout,'on',[],0,1);
colormap(cmap_ttest);
colorbar;
title('Effect of SW on RT')
caxis([-1 1]*5)
if ~isempty(find(coeff_RT(matching_elec,2)<fdr(coeff_RT(matching_elec,2),0.05)))
    ft_plot_lay_me(layout, 'chanindx',find(coeff_RT(matching_elec,2)<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
end