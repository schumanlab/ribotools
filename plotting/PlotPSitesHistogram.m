%% PlotPSiteHistogram
clc
clear variables
close all

fh = fopen('/Users/tushevg/Desktop/20190212_Npl_Mono_psiteProjection.txt','r');
txt = textscan(fh, '%s %s %s','delimiter','\t');
fclose(fh);
fname = txt{1};
val = txt{2};
cnt = txt{3};

Y = zeros(length(val), 300);
idxStart = (1:100);
idxCenter = (101:200);
idxEnd = (201:300);

Yc = zeros(length(val), 300);


for k = 1 : length(val)
    
    v = sscanf(val{k}, '%f,');
    c = sscanf(cnt{k}, '%f,');
    Y(k,idxStart) = v(idxStart)./mean(v(idxStart));
    Y(k,idxCenter) = v(idxCenter)./mean(v(idxCenter));
    Y(k,idxEnd) = v(idxEnd)./mean(v(idxEnd));
    
    Yc(k,idxStart) = c(idxStart)./mean(c(idxStart));
    Yc(k,idxCenter) = c(idxCenter)./mean(c(idxCenter));
    Yc(k,idxEnd) = c(idxEnd)./mean(c(idxEnd));
    
end


tmp = sum(Y);

tmp = tmp(26:end-25);

A = sum(tmp(1:3:end));
B = sum(tmp(2:3:end));
C = sum(tmp(3:3:end));



