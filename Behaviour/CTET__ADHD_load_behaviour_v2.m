%%
clear all;
close all;

run ../localdef_ADHD_CTET.m
% addpath(path_fieldtrip);
% ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath(path_RainCloudPlot));
%%%

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

Hit_RT_CTR=[];
for nc=1:length(ctrs)
Hit_RT_CTR(nc)=nanmean(table.Hit_RT(table.SubID==ctrs(nc)));
end
Hit_RT_ADHD=[];
for nc=1:length(adhds)
Hit_RT_ADHD(nc)=nanmean(table.Hit_RT(table.SubID==adhds(nc)));
end

stdRT_CTR=[];
for nc=1:length(ctrs)
stdRT_CTR(nc)=nanmean(table.stdRT(table.SubID==ctrs(nc)))./nanmean(table.Hit_RT(table.SubID==ctrs(nc)));
end
stdRT_ADHD=[];
for nc=1:length(adhds)
stdRT_ADHD(nc)=nanmean(table.stdRT(table.SubID==adhds(nc)))./nanmean(table.Hit_RT(table.SubID==adhds(nc)));
end

%%BarPlots
figure;
simpleBarPlot(1,100*Miss_CTR, Colors(1,:),0.9,'k',[],3);
simpleBarPlot(2,100*Miss_ADHD, Colors(2,:),0.9,'k',[],3);
format_fig; title('Miss');
set(gca,'XTick',1:2,'XTickLabel',{'CTR','ADHD'});

%%RainCloudPlots Miss

figure;
h1 = raincloud_plot(Miss_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(Miss_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'CTR', 'ADHD'});
title(['Miss']);
set(gca,'XLim', [-0.7 1.7], 'YLim', [-0.6 2]);
box off

%%RainCloudPlots False Alarms

figure;
h1 = raincloud_plot(FA_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(FA_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'CTR', 'ADHD'});
title(['False Alarms']);
set(gca,'XLim', [-0.02 0.07], 'YLim', [-40 80]);
box off

%%RainCloudPlots Hit_RT

figure;
h1 = raincloud_plot(Hit_RT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(Hit_RT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'CTR', 'ADHD'});
title(['Hit RT']);
set(gca,'XLim', [0.9 2.1], 'YLim', [-2 6]);
box off

%%RainCloudPlots std_RT

figure;
h1 = raincloud_plot(stdRT_CTR, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(stdRT_ADHD, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'CTR', 'ADHD'});
title(['stdRT']);
set(gca,'XLim', [-0.1 0.4], 'YLim', [-4 10]);
box off