%% calculateElongationRate
clc
close all

%{
clear variables
load dataTracks_TimePoints.mat;
%}
clearvars -except bed;

%% read normalization factors
fh = fopen('/Users/tushevg/Desktop/librarySize_16Apr2019.txt', 'r');
txt = textscan(fh, '%n %n %n', 'delimiter', '\t');
fclose(fh);
metainfo(:,1) = txt{1};
metainfo(:,2) = txt{2};
metainfo(:,3) = txt{3};

nFiles = length(metainfo);
nGenes = length(bed);
hFilter = ones(1, 5);
g = 1044;
%g = 10;
idxShort = metainfo(:,2) <= 3;
idxLong = metainfo(:,2) > 3;

[xShort, ~, idxRepShort] = unique(metainfo(idxShort, 1));
[xLong, ~, idxRepLong] = unique(metainfo(idxLong, 1));

ratios = ones(length(metainfo), 1);
ratios(idxShort) = mean(metainfo(idxShort, 3)) ./ metainfo(idxShort, 3);
ratios(idxLong) = mean(metainfo(idxLong, 3)) ./ metainfo(idxLong, 3);


idxThree = nchoosek(1:6, 3);
idxFour = nchoosek(1:6, 4);
idxFive = nchoosek(1:6, 5);
idxSix = nchoosek(1:6, 6);


threshX_left = 40;
threshX_right = 50;
normX = 50;

result = zeros(nGenes, 5);

for g = 1 : nGenes

    % cds coordinates
    cdsStart = bed(g).cdsStart;
    cdsEnd = bed(g).cdsEnd;
    cdsSpan = cdsEnd - cdsStart;
    offset = mod(cdsSpan, 3);
    cdsEnd = cdsEnd - offset;
    cdsSpan = cdsEnd - cdsStart;
    
    % extract tracks
    track = bed(g).linearCoverage + 1;
    track = double(track(:,cdsStart + 1 : cdsEnd));
    cdsSpanCodons = cdsSpan/3;
    
    %if cdsSpanCodons < 200
    %    continue;
    %end
    
    
    data = zeros(nFiles, cdsSpanCodons);
    ndata = zeros(nFiles, cdsSpanCodons - (threshX_left + threshX_right) + 1);
    
    xbin = (1:cdsSpanCodons);
    idxuse = (threshX_left <= xbin) & (xbin <= cdsSpanCodons - threshX_right);
    
    for f = 1 : nFiles
        trackCodons = track(f,:);
        trackCodons = sum(reshape(trackCodons, 3, length(trackCodons)/3), 1);
        trackCodons = conv(trackCodons, hFilter, 'same');
        data(f,:) = trackCodons;
        %normTrackCodons = trackCodons(idxuse);
        %normTrackCodons = cumsum(normTrackCodons ./ mean(normTrackCodons(end-normX:end)));
        %ndata(f,:) = normTrackCodons./sum(idxuse);
    end
    
    
    tmp = bsxfun(@rdivide, data, mean(data(end-50:end), 2));
    
    
    % collapse replica
    data0_short = mean(tmp(1:3,:), 1);
    data15_short = mean(tmp(7:9,:), 1);
    data30_short = mean(tmp(10:12,:), 1);
    data45_short = mean(tmp(13:15,:), 1);
    
    data0_long = mean(tmp(4:6,:), 1);
    data90_long = mean(tmp(16:18,:),1);
    data120_long = mean(tmp(19:21,:),1);
    data150_long = mean(tmp(22:24,:),1);
    
    data15_norm = data15_short ./ data0_short;
    data30_norm = data30_short ./ data0_short;
    data45_norm = data45_short ./ data0_short;
    data90_norm = data90_long ./ data0_long;
    data120_norm = data120_long ./ data0_long;
    data150_norm = data150_long ./ data0_long;
    
    %[~,x] = min(abs([data15_norm;data30_norm;data45_norm;data90_norm;data120_norm;data150_norm] - 1), [], 2);
    
    t = [15;30;45;90;120;150];
    x(1,1) = max([1,find(data15_norm <= 0.8, 1, 'last')]);
    x(2,1) = max([1,find(data30_norm <= 0.8, 1, 'last')]);
    x(3,1) = max([1,find(data45_norm <= 0.8, 1, 'last')]);
    x(4,1) = max([1,find(data90_norm <= 0.8, 1, 'last')]);
    x(5,1) = max([1,find(data120_norm <= 0.8, 1, 'last')]);
    x(6,1) = max([1,find(data150_norm <= 0.8, 1, 'last')]);
    
    
    
    [bp, cp] = bootstrapLinReg(t, x, idxThree);
    
    
    
    
    %figure('color','w');
    %plot(t,x, 'k.','markersize',18)
    
    [b, stats] = robustfit(t, x);
    result(g,1) = cdsSpanCodons;
    result(g,2) = b(2);
    result(g,3) = corr(t, x);
    result(g,4) = bp;
    result(g,5) = cp;
    
    
    %{
    figure('color','w');
    hold on;
    h(1) = plot(data15_short ./ data0_short);
    h(2) = plot(data30_short ./ data0_short);
    h(3) = plot(data45_short ./ data0_short);
    h(4) = plot(data90_long ./ data0_long);
    h(5) = plot(data120_long ./ data0_long);
    h(6) = plot(data150_long ./ data0_long);
    hold off;
    hl = legend(h,'15','30','45','90','120','150');
    %}
    
    
    %{
    % collapse replica
    data0_short = mean(ndata(1:3,:), 1);
    data15_short = mean(ndata(7:9,:), 1);
    data30_short = mean(ndata(10:12,:), 1);
    data45_short = mean(ndata(13:15,:), 1);
    thresh_short = 0.5 * data45_short(end);
    
    data0_long = mean(ndata(4:6,:), 1);
    data90_long = mean(ndata(16:18,:),1);
    data120_long = mean(ndata(19:21,:),1);
    data150_long = mean(ndata(22:24,:),1);
    thresh_long = 0.5 * data150_long(end);
    
    
    [~, yShort] = min(abs([data0_short;...
                         data15_short;...
                         data30_short;...
                         data45_short] - 0.8), [], 2);
    [~, yLong] = min(abs([data0_long;...
                         data90_long;...
                         data120_long;...
                         data150_long] - 0.8), [], 2);
    %}                
    
    %{
    figure('color','w');
    hold on;
    plot(data0_short, 'k');
    plot(data15_short, 'r');
    plot(data30_short, 'g');
    plot(data45_short, 'b');
    %plot([0,cdsSpanCodons],[thresh_short, thresh_short], 'color', [.65,.65,.65]);
    hold off;
    
    figure('color','w');
    hold on;
    plot(data0_long, 'k');
    plot(data90_long, 'r');
    plot(data120_long, 'g');
    plot(data150_long, 'b');
    %plot([0, cdsSpanCodons],[thresh_long, thresh_long], 'color', [.65,.65,.65]);
    hold off;
    %}
    
    %{
    yShortDelta = yShort - yShort(1);
    yLongDelta = yLong - yLong(1);
    tDelta = [15;30;45;90;120;150];
    yDelta = [yShortDelta(2:end);yLongDelta(2:end)];
    %}
    
    %{
    pShort = polyfit(xShort, yShort, 1);
    cShort = corr(xShort, yShort).^2;
    pLong = polyfit(xLong, yLong, 1);
    cLong = corr(xLong, yLong).^2;
    
    result(g,1) = sum(idxuse);
    result(g,2) = pShort(1);
    result(g,3) = pLong(1);
    result(g,4) = cShort;
    result(g,5) = cLong;
    %}
    

end

%{
idxsht = (result(:,1) > 400-140) & (result(:,2) > 0) & (result(:,4) >= 0.8);
idxlng = (result(:,1) > 400-140) & (result(:,3) > 0) & (result(:,5) >= 0.8);

figure('color','w');
plot(log10(result(idxsht,1)), result(idxsht, 2),'.');

figure('color','w');
plot(log10(result(idxlng,1)), result(idxlng, 3),'.');
%}

function [b,c] = bootstrapLinReg(x,y,idx)

    n = size(idx,1);
    b = zeros(n, 1);
    c = zeros(n, 1);
    for k = 1 : n
        
        xnow = x(idx(k,:));
        ynow = y(idx(k,:));
        p = polyfit(xnow, ynow, 1);
        c(k) = corr(xnow, ynow);
        b(k) = p(1);   
    end
    
    [c, i] = max(c);
    b = b(i);
end
