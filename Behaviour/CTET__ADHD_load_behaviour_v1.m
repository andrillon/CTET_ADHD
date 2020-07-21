%%
clear all;
close all;

run ../localdef_ADHD_CTET.m
% addpath(path_fieldtrip);
% ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath(path_RainCloudPlot));

%%
table=readtable([save_path filesep 'CTET_ADHD_behav_resblock.txt']);
table.SubID=categorical(table.SubID);
table.Group=categorical(table.Group);
table.Group=reordercats(table.Group,[2,1]);
table.stdRTrel=table.stdRT./table.Hit_RT;

mdl1=fitlme(table,'Hit_RT~1+Group*BlockN+(1|SubID)');
mdl2=fitlme(table,'FA~1+Group*BlockN+(1|SubID)');
mdl3=fitlme(table,'Miss~1+Group*BlockN+(1|SubID)');
mdl4=fitlme(table,'stdRT~1+Group*BlockN+(1|SubID)');
%%
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

ctrs=unique(table.SubID(table.Group=='CTR' ));
Miss_CTR=[];
for nc=1:length(ctrs)
Miss_CTR(nc)=nanmean(table.Miss(table.SubID==ctrs(nc)));
end
adhds=unique(table.SubID(table.Group=='ADHD' ));
Miss_ADHD=[];
for nc=1:length(adhds)
Miss_ADHD(nc)=nanmean(table.Miss(table.SubID==adhds(nc)));
end

FA_CTR=[];
for nc=1:length(ctrs)
FA_CTR(nc)=nanmean(table.FA(table.SubID==ctrs(nc)));
end
FA_ADHD=[];
for nc=1:length(adhds)
FA_ADHD(nc)=nanmean(table.FA(table.SubID==adhds(nc)));
end

stdRT_CTR=[];
for nc=1:length(ctrs)
stdRT_CTR(nc)=nanmean(table.stdRT(table.SubID==ctrs(nc)))./nanmean(table.Hit_RT(table.SubID==ctrs(nc)));
end
stdRT_ADHD=[];
for nc=1:length(adhds)
stdRT_ADHD(nc)=nanmean(table.stdRT(table.SubID==adhds(nc)))./nanmean(table.Hit_RT(table.SubID==adhds(nc)));
end

figure;
h1 = raincloud_plot(100*Miss_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',8);
h2 = raincloud_plot(100*Miss_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',8);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0 100], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('MISS');

figure;
h1 = raincloud_plot(100*FA_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.5);
h2 = raincloud_plot(100*FA_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.5);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0 6], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('FA');

%%
figure;
h1 = raincloud_plot(stdRT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.02);
h2 = raincloud_plot(stdRT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.02);

set(h1{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(h2{2},'LineWidth',2,'SizeData',72,'MarkerFaceAlpha',0.7);
set(gca,'XLim', [0.0 0.3], 'YLim', ylim.*[0.5 1.7]);
format_fig; title('stdRT/meanRT');
