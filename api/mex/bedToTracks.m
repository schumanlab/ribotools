%% bedToTracks
clc
close all

clear variables
bed = readBedFile('/Users/tushevg/Desktop/data/bed/ncbiRefSeq_rn6_Sep2018_appris.bed');
name = {bed.transcript}';
gene = {bed.gene}';
cdsStart = [bed.cdsStart]';
cdsSpan = [bed.cdsSpan]';
cdsSpan = cdsSpan - mod(cdsSpan, 3);

%% bed files
fileList = ls('/Users/tushevg/Desktop/data/timepoints/*.bed.gz');
fileList(end) = [];
fileList = regexp(fileList, '\n', 'split')';
nFiles = length(fileList);
fileName = regexprep(fileList, '/Users/tushevg/Desktop/data/timepoints/','');
fileName = regexprep(fileName, '.bed.gz', ''); 
fileName = cellfun(@(x) {sscanf(x,'%ds_%d')'},fileName);
metainfo = cell2mat(fileName);


%% allocate handels
handles = zeros(nFiles, 1, 'uint64');
for k = 1 : nFiles
    handles(k) = coverage('new', fileList{k});
end

%% track matrix
txName = 'NM_012920.1';
idx = strcmp(txName, name);
%k = find(idx,1,'first');
k = 86;
i = 0;
%fw = fopen('/Users/tushevg/Desktop/temp.txt','w');
tic;
for k = 1 : length(bed)
    
    qName = name{k};
    qStart = cdsStart(k);
    qEnd = cdsStart(k) + cdsSpan(k);%cdsStart(k) + cdsSpan(k); %bed(k).geneSpan;%
    
    %fprintf('%s:%d-%d\n',name{k},qStart,qEnd);
    tracks = zeros(nFiles, qEnd - qStart, 'int32');
    counts = zeros(nFiles, 1, 'int32');
    for f = 1 : nFiles
        [tracks(f,:), counts(f,1)] = coverage('query', handles(f), qName, qStart, qEnd);
    end
    
    if 15*sum(counts(1:6)) / (6 * cdsSpan(k)) < 2
        continue;
    end
    
    %
    i = i + 1;
    data(i,1).name = name{k};
    data(i,1).gene = gene{k};
    %data(i,1).cdsStart = cdsStart(k);
    %data(i,1).cdsLength = cdsSpan(k);
    data(i,1).counts = counts;
    data(i,1).tracks = tracks;
    %}
end
%fclose(fw);
toc
disp(i);

%% deallocate handels
for k = 1 : nFiles
    coverage('delete', handles(k));
end

%{

for k = 1 : length(bed)
    data(k,1).name = name{k};
    data(k,1).gene = gene{k};
    data(k,1).cdsStart = cdsStart(k);
    data(k,1).cdsLength = cdsSpan(k);
    data(k,1).tracks = zeros(cdsSpan(k), 24, 'int32');
end
save('/Users/tushevg/Desktop/test.mat','data');
%}

%{
%% fileList
fileList = ls('/Users/tushevg/Desktop/data/timepoints/*.gbed.gz');
fileList(end) = [];
fileList = regexp(fileList, '\n', 'split')';
nFiles = length(fileList);
fileName = regexprep(fileList, '/Users/tushevg/Desktop/data/timepoints/','');
fileName = regexprep(fileName, '_transcriptome_sorted_umi.gbed.gz', ''); 
fileName = cellfun(@(x) {sscanf(x,'%ds_%d')'},fileName);
metainfo = cell2mat(fileName);

%% allocate handels
handles = zeros(nFiles, 1, 'uint64');
for k = 1 : nFiles
    handles(k) = coverage('new', fileList{k});
end

%% track matrix
txName = 'NM_012920.1';
idx = strcmp(txName, name);
k = 1;
tic
for k = 1 : length(bed)
    
    tracks = zeros(nFiles, bed(k).geneSpan, 'int32');
    
    for f = 1 : nFiles
        tracks(f,:) = coverage('query', handles(f), bed(k).transcript, 0, bed(k).geneSpan);
    end
    
    bed(k).linearCoverage = tracks;

end
toc
save('dataTracks_TimePoints.mat','bed', 'metainfo')
%}


%{
txStart = 0;
txEnd = 1498;
txSpan = txEnd - txStart;

%}


%{
%% deallocate handels
for k = 1 : nFiles
    coverage('delete', handles(k));
end

%}






%% --- FUNCTIONS --- %%
function bed = readBedFile(fileBed)

    fh = fopen(fileBed, 'r');
    txt = textscan(fh, '%s', 'delimiter', '\n');
    fclose(fh);
    
    raw = txt{1};
    bed = repmat(BedLine(), length(raw), 1);
    
    for b = 1 : length(raw)
        tmp = BedLine();
        tmp.parse(raw{b});
        bed(b) = tmp;
    end

end

