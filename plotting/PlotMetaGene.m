%% PlotMetaGene
clc
clear variables
close all

%% load table
file = 'testMeta.txt';
fh = fopen(file, 'r');
hdr = regexp(fgetl(fh),'\t','split');
hdr = hdr(2:end);
fmt = repmat({'%n'}, length(hdr) + 1, 1);
fmt = sprintf('%s ', fmt{:});
fmt(end) = [];
txt = textscan(fh, fmt, 'delimiter', '\t');
fclose(fh);
xbin = txt{1};
values = [txt{2:end}];

%% plot data
idx = (4:6);
y = mean(values(:, idx), 2);
dy = std(values(:, idx), [], 2)./sqrt(length(idx));


xpatch = [xbin', fliplr(xbin')];
ypatch = [(y-dy)', fliplr((y+dy)')];
xtick = -0.2:0.2:1.2;
xticklabel = num2cell(xtick);
xticklabel(2) = {'start'};
xticklabel(end-1) = {'stop'};
figure('color','w');
hold(gca,'on');
hp = fill(xpatch, ypatch, [.75, .75, .75]);
plot(xbin, y, 'k', 'linewidth', 1.2);
set(hp, 'EdgeColor', 'none');
hold(gca,'off');
set(gca,'box','off',...
        'xlim',[-0.3,1.3],...
        'xtick', xtick,...
        'xticklabel',xticklabel,...
        'fontsize',16,...
        'ylim',[0,2]);
xlabel('relative position','fontsize',16);
ylabel('relative depth','fontsize',16);
title('Somata Total 1000','fontsize',16,'fontweight','normal');
print(gcf,'-dsvg','-r300','figureMeta_SomataTotal.svg');



