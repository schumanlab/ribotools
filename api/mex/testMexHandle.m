%% testMexHandle
clc
clear variables
close all

%% compile
mex coverage.cpp -I/usr/local/include -L/usr/local/lib -lhts;

%% execute
obj = coverage('new', '/Users/tushevg/Desktop/data/timepoints/000s_01_transcriptome_sorted_umi.gbed.gz');
x = coverage('query', obj, 200, 100);
coverage('delete', obj);