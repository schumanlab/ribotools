%% algorithmInitiationRate
clc
close all

%{
clear variables
load dataTracks_TimePoints.mat;
%}
clearvars -except bed;

%% bed names
names = {bed.gene}';

%% read CPM
fh = fopen('/Users/tushevg/Desktop/20190401_elongation_rate_total_rna_cpm.txt', 'r');
fgetl(fh);
txt = textscan(fh, '%s %n %n %n %n %n %n %n %n %n', 'delimiter', '\t');
fclose(fh);
cpmNames = txt{1};
cpmData = [txt{2:end}];
cpmData0 = mean(cpmData(:,4:6), 2);
cpmData150 = mean(cpmData(:,7:9), 2);



%% read metainfo
fh = fopen('/Users/tushevg/Desktop/librarySize_16Apr2019.txt', 'r');
txt = textscan(fh, '%n %n %n', 'delimiter', '\t');
fclose(fh);
metaTime = txt{1};
metaReplica = txt{2};
metaSizes = txt{3};
metaFactors = mean(metaSizes) ./ metaSizes;

%% meta indexes
idxShort = metaReplica <= 3;
idxLong = metaReplica > 3;
[xShort, ~, idxRepShort] = unique(metaTime(idxShort));
[xLong, ~, idxRepLong] = unique(metaTime(idxLong));

%% size
countFiles = size(metaTime, 1);
countGenes = size(bed, 1);

%% loop over genes
inirate = zeros(countGenes,2);
g = 26;
for g = 1 : countGenes

    idxCpm = strcmp(bed(g).gene, cpmNames);
    tmpCpm = [cpmData0(idxCpm);cpmData150(idxCpm)] + 1;
    cpmFactor = mean(tmpCpm) ./ tmpCpm;

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

    data = bsxfun(@times, track_codon, metaFactors);
    
    data0_short = mean(data(1:3,:));
    data0_long = mean(data(4:6,:));
    data15 = mean(data(7:9,:));
    data30 = mean(data(10:12,:));
    data45 = mean(data(13:15,:));
    data90 = mean(data(16:18,:));
    data120 = mean(data(19:21,:));
    data150 = mean(data(22:24,:));
    
    iniCodon = [data0_long(1)+1;data150(1)+1].*cpmFactor;
    inirate(g,1) = diff(iniCodon)./150;
    inirate(g,2) = span_codon;
    
end

%% read query list
fh =  fopen('/Users/tushevg/Desktop/20190302_neuropil_go_elongating_poly.txt', 'r');
txt = textscan(fh, '%s', 'delimiter', '\t');
fclose(fh);
list = txt{1};

name = {bed.gene}';
cdsSpan = [bed.cdsSpan]';

[~,idxIntersect] = intersect(name, list);
idxList = false(length(name), 1);
idxList(idxIntersect) = true;



