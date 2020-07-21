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

%%
figure;
simpleBarPlot(1,100*Miss_CTR, Colors(1,:),0.9,'k',[],3);
simpleBarPlot(2,100*Miss_ADHD, Colors(2,:),0.9,'k',[],3);
format_fig; title('MISS');
set(gca,'XTick',1:2,'XTickLabel',{'CTR','ADHD'});
