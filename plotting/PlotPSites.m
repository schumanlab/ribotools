%% PlotPSites
clc
clear variables
close all


%% read data
data = readPSiteTable('/Users/tushevg/Desktop/psiteProject_SomataMonoEnriched_01Apr2019.txt');
idx = false(length(data.file), 1);
idx(1:3) = true;


%% Plot Psites
Y = data.position(idx,:);
Y = Y ./ mean(Y(:));
Yavg = mean(Y, 1);
Yste = std(Y, [], 1) ./ sqrt(sum(idx));
X = (-25:75);

xTickLabel = arrayfun(@(x) {sprintf('%d',x)},[(-25:25:75),(-50:25:50),(-75:25:25)]);
xTickLabel(2) = {'start'};
xTickLabel(8) = {'center'};
xTickLabel(14) = {'end'};

figure('Color','w');
hold on;
xTickStart = plotTrack(X, Yavg(1:101), Yste(1:101));
xTickCenter = plotTrack(X + 122, Yavg(102:202), Yste(102:202));
xTickEnd = plotTrack(X + 244, Yavg(203:303), Yste(203:303));
hold off;
set(gca, 'box','off',...
        'xtick',[xTickStart, xTickCenter, xTickEnd],...
        'xticklabel',xTickLabel,...
        'ylim',[0,9],...
        'ytick',(0:8));
xlabel('relative position');
ylabel('normalized p-site coverage');


%% Plot Frame
pos = get(gca,'Position');
pos(1) = pos(1) + pos(3)/1.9;
pos(2) = pos(2) + pos(4)/1.6;
pos(3) = pos(3)/4;
pos(4) = pos(4)/3;


result = (data.frame(idx,:)-(1/3))./(1/3);
resultAvg = mean(result, 1);
resultSte = std(result,[],1)./sqrt(sum(idx));
p = anova1(result,[],'off');

ax = axes('Position',pos);
hold(ax,'on');
bar(ax, 0, resultAvg(1), 0.4,'facecolor',[.1,.1,.1]);
bar(ax, [1,2], resultAvg(2:3),0.4,'facecolor',[.45,.45,.45]);
plot([0,0],[resultAvg(1),resultAvg(1)+resultSte(1)],'k','linewidth',1.2);
plot([1,1],[resultAvg(2),resultAvg(2)-resultSte(2)],'k','linewidth',1.2);
plot([2,2],[resultAvg(3),resultAvg(3)-resultSte(3)],'k','linewidth',1.2);

hold(ax,'off');
set(ax,'xtick',(0:2),'xlim',[-0.5,2.5],'ylim',[-0.6,0.6],'ytick',(-0.5:0.5:0.5));
text(ax,1,0.6,sprintf('anova p = %.4e',p),'verticalalignment','bottom','horizontalalignment','center');
xlabel(ax,'frame');
ylabel(ax,'(observed - expected) / expected');


print(gcf,'-dsvg','-r300','/Users/tushevg/Desktop/figurePeriodicity_SomataMonoEnriched_01Apr2019.svg');




%% FUNCTIONS
function xTick = plotTrack(X, Yavg, Yste)
    
    Ctrack = [.5, .5, .5];
    Cframe0 = [.1, .1, .1];
    Cframe1 = [.45, .45, .45];
    Cframe2 = [.45, .45, .45];
    
    Xpatch = [X, fliplr(X)];
    Ypatch = [(Yavg - Yste), fliplr(Yavg + Yste)];
    hp = fill(Xpatch, Ypatch, Ctrack);
    set(hp, 'edgecolor', 'none', 'facealpha', .2);
    plot(X, Yavg, 'color', Ctrack);
    idxFrame0 = mod(X,3) == 0;
    idxFrame1 = mod(X,3) == 1;
    idxFrame2 = mod(X,3) == 2;
    
    plot(X(idxFrame0), Yavg(idxFrame0), '.', 'markersize', 5, 'color', Cframe0);
    plot(X(idxFrame1), Yavg(idxFrame1), '.', 'markersize', 3, 'color', Cframe1);
    plot(X(idxFrame2), Yavg(idxFrame2), '.', 'markersize', 3, 'color', Cframe2);
    xTick = (X(1):25:X(end));

end



function data = readPSiteTable(fileName)

    fh = fopen(fileName, 'r');
    txt = textscan(fh, '%s %s %s %s %s', 'delimiter', '\t');
    fclose(fh);
    data.file = txt{1};
    frame = txt{2};
    positionStart = txt{3};
    positionCenter = txt{4};
    positionEnd = txt{5};
    
    n = length(data.file);
    data.frame = zeros(n, 3);
    data.position = zeros(n, 303);
    
    
    for k = 1 : n
        data.frame(k,:) = sscanf(frame{k}, '%f,');
        data.position(k,1:101) = sscanf(positionStart{k}, '%f,');
        data.position(k,102:202) = sscanf(positionCenter{k}, '%f,');
        data.position(k,203:303) = sscanf(positionEnd{k}, '%f,');
    end
    
    %% normalize
    data.frame = bsxfun(@rdivide, data.frame, sum(data.frame, 2));
    data.position = bsxfun(@rdivide, data.position, mean(data.position, 2));
    
end