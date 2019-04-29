%% algorithmElongationRate
clc
close all

%{
clear variables
load dataTracks_TimePoints.mat;
%}
clearvars -except bed;

%% read query list
fh =  fopen('/Users/tushevg/Desktop/20190302_neuropil_go_elongating_mono.txt', 'r');
txt = textscan(fh, '%s', 'delimiter', '\t');
fclose(fh);
list = txt{1};

name = {bed.gene}';
cdsSpan = [bed.cdsSpan]';

[~,idxIntersect] = intersect(name, list);
idxList = false(length(name), 1);
idxList(idxIntersect) = true;



%% read metainfo
fh = fopen('/Users/tushevg/Desktop/librarySize_16Apr2019.txt', 'r');
txt = textscan(fh, '%n %n %n', 'delimiter', '\t');
fclose(fh);
metaTime = txt{1};
metaReplica = txt{2};
metaSizes = txt{3};

%% meta indexes
idxShort = metaReplica <= 3;
idxLong = metaReplica > 3;
[xShort, ~, idxRepShort] = unique(metaTime(idxShort));
[xLong, ~, idxRepLong] = unique(metaTime(idxLong));

%% size
countFiles = size(metaTime, 1);
countGenes = size(bed, 1);

%% loop over genes
dataA = zeros(3, 650, countGenes);
countsA = false(countGenes,1);
dataB = zeros(3, 650, countGenes);
countsB = false(countGenes,1);
g = 1;
for g = 1 : countGenes

    % extract nucleotide track
    offset = mod(bed(g).cdsSpan, 3);
    track_nucleotide = double(bed(g).linearCoverage(:, (bed(g).cdsStart + 1):(bed(g).cdsEnd - offset)) + 1);
    
    % extract codon track
    span_codon = size(track_nucleotide, 2) / 3;
    track_codon = zeros(countFiles, span_codon);
    for f = 1 : countFiles
        track_codon(f,:) = sum(reshape(track_nucleotide(f,:), 3, span_codon), 1);
        track_codon(f,:) = movmean(track_codon(f,:), 5);
    end

    data = track_codon;
    
    data0_short = mean(data(1:3,:));
    data0_long = mean(data(4:6,:));
    data15 = mean(data(7:9,:));
    data30 = mean(data(10:12,:));
    data45 = mean(data(13:15,:));
    data90 = mean(data(16:18,:));
    data120 = mean(data(19:21,:));
    data150 = mean(data(22:24,:));
    
    mtxA = [data15./data0_short;...
            data30./data0_short;...
            data45./data0_short];
   
    mtxB = [data45./data0_short;...
            data90./data0_long;...
            data120./data0_long];
    
    
    %% split by length
    if (400 <= span_codon)
        
        nrmAa = bsxfun(@rdivide, mtxA, mean(mtxA(:,350:400), 2));
        %nrmAb = bsxfun(@rdivide, mtxA, mean(mtxA(:,end-50:end), 2));
        %nrmAc = bsxfun(@rdivide, mtxA, mean(mtxA(:,end-75:end-25), 2));
        
        j = min([550,size(nrmAa,2)]);
        dataA(:,1:j,g) = dataA(:,1:j,g) + nrmAa(:,1:j);
        %dataA(:,:,2) = dataA(:,:,2) + nrmAb(:,1:400);
        %dataA(:,:,3) = dataA(:,:,3) + nrmAc(:,1:400);
        
        countsA(g) = true;
        
    end
    
    
    
    if (650 <= span_codon)
        
        nrmBa = bsxfun(@rdivide, mtxB, mean(mtxB(:,600:650), 2));
        %nrmBb = bsxfun(@rdivide, mtxB, mean(mtxB(:,end-50:end), 2));
        %nrmBc = bsxfun(@rdivide, mtxB, mean(mtxB(:,end-75:end-25), 2));
        
        dataB(:,:,g) = dataB(:,:,1) + nrmBa(:,1:650);
        %dataB(:,:,2) = dataB(:,:,2) + nrmBb(:,1:550);
        %dataB(:,:,3) = dataB(:,:,3) + nrmBc(:,1:550);
        
        countsB(g) = true;
        
    end
    
end

dataA = mean(dataA(:,:,countsA & idxList), 3);
dataB = mean(dataB(:,:,countsB & idxList), 3);
plotTimePoints(dataA, 1, {'15','30','45'},[15;30;45],'A', 0.9);
plotTimePoints(dataB, 1, {'45','90','120'},[45;90;120],'B', 0.8);

function plotTimePoints(data, counts, labels, t, name, thresh)
    h = zeros(size(data,1),1);
    xbin = (1:size(data,2));
    x = zeros(size(data,1), 1);
    
    figure('color','w');
    hold on;
    plot([0,size(data,2)],[1,1],'k');
    for k = 1 : size(data,1)
        h(k) = plot(data(k,:)./counts);
        x(k) = find(((data(k,:)./counts) >= thresh) & (xbin > 40), 1, 'first');
    end
    hold off;
    hl = legend(h, labels);
    set(hl,'edgecolor','w');
    set(gca,'fontsize',14,'xlim',[0,550],'ylim',[0.2,1.2]);
    xlabel('distance [codons]','fontsize',14);
    ylabel('normalized coverage','fontsize',14);
    title(name,'fontsize',14,'fontweight','normal');
    
    p = polyfit(t, x, 1);
    xfit = linspace(min(t),max(t), 100);
    yfit = polyval(p, xfit);
    figure('color','w');
    hold on;
    plot(xfit, yfit, 'color', [.45,.45,.45]);
    plot(t, x, 'k.','markersize',14);
    hold off;
    title(sprintf('%s : y(t) = %.4f + %.4f',name, p(1), p(2)),'fontsize',14,'fontweight','normal');
    set(gca,'fontsize',14,'xtick',t);
    xlabel('time [sec]','fontsize',14);
    ylabel('position [codon]','fontsize',14);
    
    
end
