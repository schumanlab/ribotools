%% alternativeElongationRate
clc
close all

%% load tracks
%{
clear variables
load /Users/tushevg/Desktop/data/timepoints/dataTracks_17Apr2019.mat
%}
clearvars -except data metainfo

timeList = metainfo(:,1);
timeList(4:6) = 1;
[timeListUnique,~,idxMerge] = unique(timeList);


timeList = metainfo(:,1);
timeList(4:6) = 1;
[timeListUnique,~,idxMerge] = unique(timeList);


nTimes = max(idxMerge);
nGenes = length(data);
g = 2; % 8359 longest


%for g = 1 : nGenes

    tracks = double(data(g).tracks);
    %counts = double(data(g).counts);
    %factor = geomean(counts)./counts;
    %tracks = bsxfun(@times, tracks, factor);
    nCodons = size(tracks, 2) / 3;
    codons = zeros(nTimes, nCodons);
    for t = 1 : nTimes
        arrayNucl = sum(tracks(idxMerge == t,:) + 1, 1);
        arrayCodons = sum(reshape(arrayNucl, 3, nCodons), 1);
        arrayCodons = movmean(arrayCodons, 5);
        arrayCodons = arrayCodons ./ mean(arrayCodons(end-ceil(0.1*nCodons):end));
        codons(t,:) = arrayCodons;
    end
    
    %% signal squered error
    sse = [(codons(3,:) - codons(1,:));...
           (codons(4,:) - codons(1,:));...
           (codons(5,:) - codons(1,:));...
           (codons(6,:) - codons(2,:));...
    	   (codons(7,:) - codons(2,:));...
           (codons(8,:) - codons(2,:))].^2;
       
    
    avg = mean(sse, 2);
    
    figure('color','w');
    plot(sse(1,:),'k');
    hold on;
    plot([0,nCodons],[avg(1),avg(1)],'r');
    hold off;
       
       
%end