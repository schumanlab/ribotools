%% testMexHandle
clc
clear variables
close all

%% compile
mex coverage.cpp -I/usr/local/include -L/usr/local/lib -lhts;

%% fileList
fileList = ls('/Users/tushevg/Desktop/*.bed.gz');
fileList(end) = [];
fileList = regexp(fileList, '\n', 'split')';
nFiles = length(fileList);

%% allocate handels
handles = zeros(nFiles, 1, 'uint64');
for k = 1 : nFiles
    handles(k) = coverage('new', fileList{k});
end

%% track matrix
txName = 'NM_012920.1';
txStart = 0;
txEnd = 1498;
txSpan = txEnd - txStart;
tracks = zeros(nFiles, txSpan, 'int32');
tic
for k = 1 : nFiles
    [tracks(k,:),tmp] = coverage('query', handles(k), txName, txStart, txEnd);
end
toc


%% deallocate handels
for k = 1 : nFiles
    coverage('delete', handles(k));
end
