%% calculateElongationRateByList
clc
close all

%% read in data
clear variables
load dataTracks_TimePoints.mat;
clearvars -except bed metainfo;

%% read in list
fh = fopen('/Users/tushevg/Desktop/20190302_neuropil_go_elongating_mono.txt', 'r');
txt = textscan(fh, '%s', 'delimiter', '\t');
fclose(fh);





