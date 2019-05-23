%% exploreElongationRate
clc
close all

%% load tracks
%{
clear variables
load /Users/tushevg/Desktop/data/timepoints/dataTracks_17Apr2019.mat
%}
clearvars -except data metainfo

listMono = readGeneList('/Users/tushevg/Desktop/data/timepoints/GeneLists/GeneListSomataMonoEnriched_20May.txt');
listPoly = readGeneList('/Users/tushevg/Desktop/data/timepoints/GeneLists/GeneListSomataPolyEnriched_20May.txt');


timeList = metainfo(:,1);
timeList(4:6) = 1;
[timeListUnique,~,idxMerge] = unique(timeList);


nTimes = max(idxMerge);
nGenes = length(data);
g = 1; % 8359 longest


arrayTime = [15;30;45;90;120;150];

nBin = 5;
X = arrayTime(1:nBin);
Xfit = [ones(nBin,1), X];


result = zeros(nTimes, 200);
resultMono = zeros(nTimes, 200);
resultPoly = zeros(nTimes, 200);
l = 0;
m = 0;
p = 0;

tic
for g = 1 : nGenes

    tracks = double(data(g).tracks);
    counts = double(data(g).counts);
    factor = geomean(counts)./counts;
    tracks = bsxfun(@times, tracks, factor);
    nCodons = size(tracks, 2) / 3;
    codons = zeros(nTimes, nCodons);
    for t = 1 : nTimes
        arrayNucl = sum(tracks(idxMerge == t,:) + 1, 1);
        arrayCodons = sum(reshape(arrayNucl, 3, nCodons), 1);
        %arrayCodons = movmean(arrayCodons, 5);
        %arrayCodons = arrayCodons ./ mean(arrayCodons(end-ceil(0.1*nCodons):end));
        codons(t,:) = arrayCodons;
    end
    
    %% signal squered error
    %{
    sse = [(codons(3,:) - codons(1,:));...
           (codons(4,:) - codons(1,:));...
           (codons(5,:) - codons(1,:));...
           (codons(6,:) - codons(2,:));...
    	   (codons(7,:) - codons(2,:));...
           (codons(8,:) - codons(2,:))];
    %}
    
    if (nCodons >= 200) && all(~isnan(codons(:)))
        
        result = result + codons(:,1:200);
        l = l + 1;
        
        if any(strcmp(data(g).gene, listMono))
            resultMono = resultMono + codons(:,1:200);
            m = m + 1;
        end
        
        if any(strcmp(data(g).gene, listPoly))
            resultPoly = resultPoly + codons(:,1:200);
            p = p + 1;
        end
        
        
    end
    
    
       
end
toc



result = bsxfun(@rdivide, cumsum(result, 2), sum(result, 2));
resultMono = bsxfun(@rdivide, cumsum(resultMono, 2), sum(resultMono, 2));
resultPoly = bsxfun(@rdivide, cumsum(resultPoly, 2), sum(resultPoly, 2));

%result = abs(bsxfun(@rdivide, cumsum(result, 2), sum(result, 2)) - 0.5);
%resultMono = abs(bsxfun(@rdivide, cumsum(resultMono, 2), sum(resultMono, 2)) - 0.5);
%resultPoly = abs(bsxfun(@rdivide, cumsum(resultPoly, 2), sum(resultPoly, 2)) - 0.5);


figure('color','w');
hold on;
plot(result(1,:),'k');
plot(result(3,:),'r');
plot(result(4,:),'g');
plot(result(5,:),'b');
hold off;

figure('color','w');
hold on;
result = resultMono;
plot(result(1,:),'k');
plot(result(3,:),'r');
plot(result(4,:),'g');
plot(result(5,:),'b');
hold off;
xlabel('polition along ORF [codons]');
ylabel('cumulative fraction of footprints');
title('Somata - Monosome');
print(gcf,'-dpng','-r300','figure_SomataMonosomeEnriched.png');

figure('color','w');
hold on;
result = resultPoly;
plot(result(1,:),'k');
plot(result(3,:),'r');
plot(result(4,:),'g');
plot(result(5,:),'b');
hold off;
xlabel('polition along ORF [codons]');
ylabel('cumulative fraction of footprints');
title('Somata - Polysome');
print(gcf,'-dpng','-r300','figure_SomataPolysomeEnriched.png');
%}


function list = readGeneList(fileName)

    fh = fopen(fileName, 'r');
    txt = textscan(fh, '%s','delimiter','\n');
    fclose(fh);
    list = txt{1};

end
