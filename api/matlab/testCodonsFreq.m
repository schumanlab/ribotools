%% testCodonsFreq
clc
close all
clear variables

fh = fopen('codonNFC_Report_28May2019.txt', 'r');
header = fgetl(fh);
txt = textscan(fh, '%s %s %s %s %n %n %n %n %n %n','delimiter','\t');
fclose(fh);
code = txt{1};
letter = txt{2};
abbreviation = txt{3};
name = txt{4};
mfdr = [txt{5:end}];
