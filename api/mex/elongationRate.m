%% elongationRate
clc
close all

%% load tracks
%{
clear variables
load /Users/tushevg/Desktop/dataTracks_17Apr2019.mat
%}
clearvars -except data metainfo

nTracks = size(data, 1);
k = 1;
i = 0;
%for k = 1 : nTracks
    
    track0_short = sum(data(k).tracks(1:3,:));
    track0_long = sum(data(k).tracks(4:6,:));
    track15 = sum(data(k).tracks(13:15,:));
    track30 = sum(data(k).tracks(16:18,:));
    track45 = sum(data(k).tracks(19:21,:));
    track90 = sum(data(k).tracks(22:24,:));
    track120 = sum(data(k).tracks(7:9,:));
    track150 = sum(data(k).tracks(10:12,:));
    
    counts = [sum(data(k).counts(1:3));...
              sum(data(k).counts(4:6));...
              sum(data(k).counts(13:15));...
              sum(data(k).counts(16:18));...
              sum(data(k).counts(19:21));...
              sum(data(k).counts(22:24));...
              sum(data(k).counts(7:9));...
              sum(data(k).counts(10:12))];
              
    fctr = mean(counts) ./ counts;
    
    tracks = [track0_short;...
              track0_long;...
              track15;...
              track30;...
              track45;...
              track90;...
              track120;...
              track150] .* fctr;
    
    nCodons = size(tracks, 2) / 3;
    if (nCodons < 400)
        %continue;
    end
    i = i + 1;
    
    covg = zeros(8, nCodons);
    for t = 1 : 8
        Yn = tracks(t,:);
        Yc = sum(reshape(Yn, 3, nCodons, 1));
        Yc = movmean(Yc, 30);
        %Yc = Yc ./ mean(Yc(end-70:end-20));
        covg(t,:) = Yc;
        
    end
    Xc = (1:nCodons)';
    
    
    figure('color','w');
    hold on;
    plot(cumsum(covg(1,:)),'k');
    plot(cumsum(covg(3,:)),'r');
    plot(cumsum(covg(4,:)),'g');
    plot(cumsum(covg(5,:)),'b');
    hold off;
    
    figure('color','w');
    hold on;
    plot(cumsum(covg(2,:)),'k');
    plot(cumsum(covg(6,:)),'r');
    plot(cumsum(covg(7,:)),'g');
    plot(cumsum(covg(8,:)),'b');
    hold off;
    
    

%end
