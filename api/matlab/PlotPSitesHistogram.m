%% PlotPSitesHistogram

clc
clear variables
close all


fr = fopen('testFrame.txt','r');
txt = textscan(fr, '%s %s %n %n %n %n %n','delimiter','\t');
fclose(fr);

files = txt{1};
labels = txt{2};
spans = txt{3};
offsets = txt{4};
counts = [txt{5:end}];
pcounts = bsxfun(@rdivide, counts, sum(counts, 2)); 

pcounts(end,:) = [];
offsets(end) = [];
spans(end) = [];
labels(end) = [];


idxLabel = strcmp('CDS',labels);
idxSpans = (29 <= spans) & (spans <= 33);
idx = idxLabel & idxSpans;

%
figure('color','w');
imagesc(pcounts(idx,:));
set(gca,'ytick',(1:sum(idx)),...
        'yticklabel',(spans(idx)),...
        'xtick',(1:3),...
        'xticklabel',(0:2));
%}