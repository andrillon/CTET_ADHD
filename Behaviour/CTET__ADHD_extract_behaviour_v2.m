%%
clear all
close all

run ../localdef_ADHD_CTET.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath(path_RainCloudPlot));


% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*' filesep '*' filesep '*_CTET*.bdf']);

wrong_events={'MC_control_CTET.bdf'};
%%
res_mat=[];
group_cond=[];
problematic_files=[];
resblock_mat=[];
groupblock_cond=[];
nc=0;
for nF=1:length(files)
    File_Name=files(nF).name;
    if ismember(File_Name,wrong_events)
        fprintf('... skipping %s\n',File_Name);
        continue;
    else
        fprintf('... processing %s\n',File_Name);
    end
    septag=findstr(File_Name,'_');
    SubN=(File_Name(1:septag(1)-1));
    if ~isempty(findstr(files(nF).folder,'adhds'))
        Group='ADHD';
    elseif ~isempty(findstr(files(nF).folder,'controls'))
        Group='CTR';
    end
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    try
    hdr=ft_read_header([files(nF).folder filesep File_Name]);
    events=ft_read_event([files(nF).folder filesep File_Name]);
    catch
                     nc=nc+1;
  problematic_files{nc}= File_Name;
               continue;
    end
    my_events=events(find_trials({events.type},'STATUS'));
    unique_values=unique([my_events.value]);
    fprintf('... %g events found:',length([my_events.value]));
    for nU=1:length(unique_values)
        fprintf(' %g (n=%g) ',unique_values(nU),sum([my_events.value]==unique_values(nU)))
    end
    fprintf('\n');


    %%% find blocks
    all_idx=[my_events.sample];
    all_val=[my_events.value];
    if isempty(all_idx)
%         pause;
        nc=nc+1;
       problematic_files{nc}= File_Name;
       continue;
    end
    blocks_boundaries=all_idx(find(diff(all_idx/hdr.Fs)>1.7)+1);
        blocks_boundaries=[all_idx(1) blocks_boundaries];
    fprintf('... %g block transitions found\n',length(blocks_boundaries));
    if length(blocks_boundaries)~=8
%         pause;
        nc=nc+1;
       problematic_files{nc}= File_Name;
       continue;
    end
    all_val(ismember(all_idx,blocks_boundaries))=100;
    
    %%% find relevant events
    resp_idx=[my_events(all_val==1).sample];
    
    targets_idx=[my_events(all_val==max(unique_values)).sample];
    nontargets_idx=[my_events(all_val==20).sample];
    stim_idx=[my_events(all_val==max(unique_values) | all_val==20).sample];
    stim_val=[my_events(all_val==max(unique_values) | all_val==20).value];
    
    %%% Gather info per targets
    this_behav=nan(length(stim_idx),9);
    for nSt=1:length(all_idx)
        this_FA=NaN;
        this_RT=NaN;
        if all_val(nSt)==max(unique_values) % we have a target
            count=0; resprec=0;
            while nSt+count<length(all_val) && all_val(nSt+count)~=1
                count=count+1;
            end
            if all_val(nSt+count)==1
                this_RT=(all_idx(nSt+count)-all_idx(nSt))/hdr.Fs;
            else
                this_RT=NaN;
            end
        elseif all_val(nSt)==20 % we have a non-target
            if nSt<length(all_val) & nSt>5
                if all_val(nSt+1)==1 && max(all_val(nSt+(-5:-1)))~=max(unique_values)
                    this_FA=1;
                else
                    this_FA=0;
                end
                count=0; resprec=0;
                while nSt+count<length(all_val) && all_val(nSt+count)~=1
                    count=count+1;
                end
                if all_val(nSt+count)==1
                    this_RT=(all_idx(nSt+count)-all_idx(nSt))/hdr.Fs;
                else
                    this_RT=NaN;
                end
            end
        else
            continue;
        end
        %%%% sort things out
        this_block=max(find(all_idx(nSt)>blocks_boundaries));
        if isempty(this_block)
            this_block=NaN;
        end
        this_behav(find(stim_idx==all_idx(nSt)),:)=[find(stim_idx==all_idx(nSt)) stim_idx(find(stim_idx==all_idx(nSt))) this_block all_val(nSt)==max(unique_values) this_RT NaN NaN NaN NaN];
    end
    this_behav(:,8)=[diff(this_behav(:,2))/hdr.Fs ; NaN];
    this_behav(this_behav(:,5)>this_behav(:,8),5)=NaN;
    
    respidx=find(~isnan(this_behav(:,5)));
    this_behav(this_behav(:,4)==0,6)=1;
    this_behav(this_behav(:,4)==1,7)=0;
    for m=1:length(respidx)
        if this_behav(respidx(m),4)==1 % resp on target
            if this_behav(respidx(m),5)>min(this_behav(:,8))
                this_behav(respidx(m),7)=1;
                this_behav(respidx(m),9)=this_behav(respidx(m),5);
            else
                this_behav(respidx(m),7)=0;
                this_behav(respidx(m)-1,6)=0;
                this_behav(respidx(m)-1,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8);
            end
        elseif this_behav(respidx(m),4)==0 % resp on non target
            if this_behav(respidx(m)-1,4)==1
                this_behav(respidx(m)-1,7)=1;
                this_behav(respidx(m)-1,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8);
            else
                if this_behav(respidx(m)-2,4)==1
                this_behav(respidx(m)-2,7)=1;
                this_behav(respidx(m)-2,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8)+this_behav(respidx(m)-2,8);    
                else
                    this_behav(respidx(m)-1,6)=0;
                    this_behav(respidx(m)-1,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8);
                end
            end
        end
    end
    
    this_cleanRT=this_behav(:,8);
    res_mat=[res_mat; [nF*ones(size(this_behav,1),1) this_behav]];
    group_cond=[group_cond ; repmat({Group},size(this_behav,1),1)];

    this_table=array2table(this_behav,...
        'VariableNames',{'TrialN','Sample','BlockN','StimType','rawRT','corrNT','corrTG','ITI','RT'});
    writetable(this_table,[save_path filesep 'CTET_ADHD_behav_' File_Name(1:end-4) '.txt']);

    for nbl=1:8
        [dprime,crit]=calc_dprime(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==1,this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==0);
        resblock_mat=[resblock_mat; [nF nbl nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==1) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==0) ...
            nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==1) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==0) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1 & this_behav(:,7)==1,9)) nanstd(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1 & this_behav(:,7)==1,9)) ...
            dprime crit]];
        groupblock_cond=[groupblock_cond ; {Group}];
    end
    
end

%% 
table1=array2table(res_mat,'VariableNames',{'SubID','TrialN','Sample','BlockN','StimType','rawRT','corrNT','corrTG','ITI','RT'});
table1.Group=group_cond;
table1.SubID=categorical(table1.SubID);
table1.Group=categorical(table1.Group);
table1.Group=reordercats(table1.Group,[2,1]);

mdl0c=fitlme(table1,'corrTG~1+BlockN+(1|SubID)');
mdl1c=fitlme(table1,'corrTG~1+BlockN+Group+(1|SubID)');
compare(mdl0c,mdl1c) % LRStat (Chi-square Likelihood Ratio Test) pValue

mdl0b=fitlme(table1,'corrNT~1+BlockN+(1|SubID)');
mdl1b=fitlme(table1,'corrNT~1+BlockN+Group+(1|SubID)');
compare(mdl0b,mdl1b) % LRStat (Chi-square Likelihood Ratio Test) pValue

mdl0a=fitlme(table1,'RT~1+BlockN+(1|SubID)');
mdl1a=fitlme(table1,'RT~1+BlockN+Group+(1|SubID)');
compare(mdl0a,mdl1a) % LRStat (Chi-square Likelihood Ratio Test) pValue


%%
table2=array2table(resblock_mat,'VariableNames',{'SubID','BlockN','CR','FA','Hit','Miss','Hit_RT','stdRT','dprime','crit'});
table2.Group=groupblock_cond;
table2.SubID=categorical(table2.SubID);
table2.Group=categorical(table2.Group);
table2.Group=reordercats(table2.Group,[2,1]);

mdl0c=fitlme(table2,'Miss~1+BlockN+(1|SubID)');
mdl1c=fitlme(table2,'Miss~1+BlockN+Group+(1|SubID)');
compare(mdl0c,mdl1c) % LRStat (Chi-square Likelihood Ratio Test) pValue

mdl0b=fitlme(table2,'FA~1+BlockN+(1|SubID)');
mdl1b=fitlme(table2,'FA~1+BlockN+Group+(1|SubID)');
compare(mdl0b,mdl1b) % LRStat (Chi-square Likelihood Ratio Test) pValue

mdl0a=fitlme(table2,'Hit_RT~1+BlockN+(1|SubID)');
mdl1a=fitlme(table2,'Hit_RT~1+BlockN+Group+(1|SubID)');
compare(mdl0a,mdl1a) % LRStat (Chi-square Likelihood Ratio Test) pValue

mdl0d=fitlme(table2,'stdRT~1+BlockN+(1|SubID)');
mdl1d=fitlme(table2,'stdRT~1+BlockN+Group+(1|SubID)');
compare(mdl0d,mdl1d) % LRStat (Chi-square Likelihood Ratio Test) pValue

% mdl0e=fitlme(table2,'Hit~1+BlockN+(1|SubID)');
% mdl1e=fitlme(table2,'Hit~1+BlockN+Group+(1|SubID)');
% compare(mdl0e,mdl1e) % LRStat (Chi-square Likelihood Ratio Test) pValue
% 
% mdl0f=fitlme(table2,'dprime~1+BlockN+(1|SubID)');
% mdl1f=fitlme(table2,'dprime~1+BlockN+Group+(1|SubID)');
% compare(mdl0f,mdl1f) % LRStat (Chi-square Likelihood Ratio Test) pValue
% 
% mdl0g=fitlme(table2,'crit~1+BlockN+(1|SubID)');
% mdl1g=fitlme(table2,'crit~1+BlockN+Group+(1|SubID)');
% compare(mdl0g,mdl1g) % LRStat (Chi-square Likelihood Ratio Test) pValue


writetable(table2,[save_path filesep 'CTET_ADHD_behav_resblock.txt']);

%%

%%
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

ctrs=unique(table2.SubID(table2.Group=='CTR' ));
Miss_CTR=[];
for nc=1:length(ctrs)
Miss_CTR(nc)=nanmean(table2.Miss(table2.SubID==ctrs(nc)));
end
adhds=unique(table2.SubID(table2.Group=='ADHD' ));
Miss_ADHD=[];
for nc=1:length(adhds)
Miss_ADHD(nc)=nanmean(table2.Miss(table2.SubID==adhds(nc)));
end

FA_CTR=[];
for nc=1:length(ctrs)
FA_CTR(nc)=nanmean(table2.FA(table2.SubID==ctrs(nc)));
end
FA_ADHD=[];
for nc=1:length(adhds)
FA_ADHD(nc)=nanmean(table2.FA(table2.SubID==adhds(nc)));
end

stdRT_CTR=[];
for nc=1:length(ctrs)
stdRT_CTR(nc)=nanmean(table2.stdRT(table2.SubID==ctrs(nc)))./nanmean(table2.Hit_RT(table2.SubID==ctrs(nc)));
end
stdRT_ADHD=[];
for nc=1:length(adhds)
stdRT_ADHD(nc)=nanmean(table2.stdRT(table2.SubID==adhds(nc)))./nanmean(table2.Hit_RT(table2.SubID==adhds(nc)));
end

Hit_RT_CTR=[];
for nc=1:length(ctrs)
Hit_RT_CTR(nc)=nanmean(table2.Hit_RT(table2.SubID==ctrs(nc)));
end
Hit_RT_ADHD=[];
for nc=1:length(adhds)
Hit_RT_ADHD(nc)=nanmean(table2.Hit_RT(table2.SubID==adhds(nc)));
end

Hit_CTR=[];
for nc=1:length(ctrs)
Hit_CTR(nc)=nanmean(table2.Hit(table2.SubID==ctrs(nc)));
end
Hit_ADHD=[];
for nc=1:length(adhds)
Hit_ADHD(nc)=nanmean(table2.Hit(table2.SubID==adhds(nc)));
end

%Miss
figure;
h1 = raincloud_plot(100*Miss_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',8,'bound_data',[0 100]);
h2 = raincloud_plot(100*Miss_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',8,'bound_data',[0 100]);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0 100], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('MISS'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

%FA
figure;
h1 = raincloud_plot(100*FA_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.5,'bound_data',[0 100]);
h2 = raincloud_plot(100*FA_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.5,'bound_data',[0 100]);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0 6], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('FA'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

%%
%stdRT
figure;
h1 = raincloud_plot(stdRT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.02);
h2 = raincloud_plot(stdRT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.02);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0.0 0.3], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('stdRT/meanRT'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

%Hit_RT
figure;
h1 = raincloud_plot(Hit_RT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.04);
h2 = raincloud_plot(Hit_RT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.04);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [1.1 1.9], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('RT (s)'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

%Hit
figure;
h1 = raincloud_plot(Hit_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.02);
h2 = raincloud_plot(Hit_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.02);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0.0 0.3], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('Hit'); legend([h1{1} h2{1}], {'Controls', 'ADHDs'});

%%
%% Repeated Measures plot
%Hit
data_to_plot=[];
group_labels={'CTR','ADHD'};
for i = 1:8 % number of repetitions
    for j = 1:2 % number of group
        data_to_plot{i, j} = table2.Hit(table2.BlockN==i & table2.Group==group_labels{j});
    end
end

figure; hold on;
h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'YLim', [-0.3 1.6]);
title(['Hits per block']);
format_fig;

%Miss
data_to_plot=[];
group_labels={'CTR','ADHD'};
for i = 1:8 % number of repetitions
    for j = 1:2 % number of group
        data_to_plot{i, j} = table2.Miss(table2.BlockN==i & table2.Group==group_labels{j});
    end
end

figure; hold on;
h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'YLim', [-0.3 1.6]);
title(['Misses per block']);
format_fig;

%False Alarms
data_to_plot=[];
group_labels={'CTR','ADHD'};
for i = 1:8 % number of repetitions
    for j = 1:2 % number of group
        data_to_plot{i, j} = table2.FA(table2.BlockN==i & table2.Group==group_labels{j});
    end
end

figure; hold on;
h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'YLim', [-0.3 1.6]);
title(['False alarms per block']);
format_fig;

%Hit RT
data_to_plot=[];
group_labels={'CTR','ADHD'};
for i = 1:8 % number of repetitions
    for j = 1:2 % number of group
        data_to_plot{i, j} = table2.Hit_RT(table2.BlockN==i & table2.Group==group_labels{j});
    end
end

figure; hold on;
h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'YLim', [-0.3 1.6]);
title(['Hit RTs per block']);
format_fig;

%stdRT
data_to_plot=[];
group_labels={'CTR','ADHD'};
for i = 1:8 % number of repetitions
    for j = 1:2 % number of group
        data_to_plot{i, j} = table2.stdRT(table2.BlockN==i & table2.Group==group_labels{j});
    end
end

figure; hold on;
h   = rm_raincloud(data_to_plot, Colors(1:2,:));
% set(gca, 'YLim', [-0.3 1.6]);
title(['stdRTs per block']);
format_fig;