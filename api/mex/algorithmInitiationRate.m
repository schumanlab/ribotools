%% algorithmInitiationRate
clc
close all

%% load tracks
%{
clear variables
load /Users/tushevg/Desktop/data/timepoints/dataTracks_17Apr2019.mat;
%}
clearvars -except data metainfo;

%% load library sizes
libsize = loadLibrarySizes('/Users/tushevg/Desktop/data/timepoints/librarySizes_03May2019.txt');

%% load RNA expression
rna = loadRNAexpression('/Users/tushevg/Desktop/data/timepoints/elongationRate_RNA_apprisCountsMatrix.txt');


%% compile matrix
%
nGenes = length(data);

counts = zeros(nGenes, 6);
codons = zeros(nGenes, 6);
avg = zeros(nGenes, 6);
rnas = zeros(nGenes,6);
L = zeros(nGenes, 1);
for g = 1 : nGenes
    name = data(g).gene;
    counts(g,:) = double([data(g).counts(4:6)',data(g).counts(10:12)']);
    tracks = double([data(g).tracks(4:6,:);
                     data(g).tracks(10:12,:)]);
    
    for t = 1 : size(tracks,1)
        
        tmp = sum(reshape(tracks(t,:), 3, size(tracks,2)/3),1);
        tmp = movmean(tmp, 5);
        codons(g,t) = sum(tmp(1:5));
        avg(g,t) = mean(tmp) + 1;
    end
    
    L(g) = size(tracks, 2)/3;
    %idxCpm = strcmp(name, rna.name);
    %rnas(g,:) = [rna.counts(idxCpm,1:3),rna.counts(idxCpm,10:12)] + 1;
end

[ncounts, factorDepth] = DESeqNormalization(counts);
ncodons = bsxfun(@rdivide, codons, factorDepth);
%factorRNA = bsxfun(@rdivide, mean(rnas, 2), rnas);
%ncodons = ncodons .* factorRNA;

nrate = codons;

%ncodons = codons ./ avg;
initRate = (ncodons(:,4:6) - ncodons(:,1:3))./150;
initRateTest = (nrate(:,4:6) - nrate(:,1:3))./150;
%initRate = max(initRate,[],2);

names = {data.gene}';
idx = strcmp('Camk2a', names);
%}
%
fw = fopen('initiationTest_DepthNorm_06May2019.txt','w');
tmp = [{data.gene}',num2cell(initRate)]';
fprintf(fw,'%s\t%.4f\t%.4f\t%.4f\n',tmp{:});
fclose(fw);
%}



%% FUNCTIONS

function cpm = loadRNAexpression(filename)

    fh = fopen(filename, 'r');
    hdr = fgetl(fh);
    hdr = regexp(hdr, '\t', 'split');
    hdr = cell2mat(cellfun(@(x) {sscanf(x,'T%ds_%d')'}, hdr(2:end))');
    txt = textscan(fh, '%s %n %n %n %n %n %n %n %n %n  %n %n %n', 'delimiter', '\t');
    fclose(fh);
    cpm.name = txt{1};
    cpm.header = hdr;
    cpm.counts = [txt{2:end}];
    [~,idxsort] = sortrows(cpm.header, [1,2]);
    cpm.counts = cpm.counts(:,idxsort);
    idxfilter = all(cpm.counts == 0, 2);
    cpm.header = cpm.header(idxsort,:);
    cpm.name(idxfilter) = [];
    cpm.counts(idxfilter,:) = [];
    
end


function libsize = loadLibrarySizes(filename)

    fh = fopen(filename, 'r');
    txt = textscan(fh, '%n %n %n', 'delimiter', '\t');
    fclose(fh);
    libsize.timepoint = txt{1};
    libsize.replica = txt{2};
    libsize.counts = txt{3};
    libsize.factor = mean(libsize.counts) ./ libsize.counts;

end

%{
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
%}


