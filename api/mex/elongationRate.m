%% elongationRate
clc
close all

%% load tracks
%{
clear variables
load /Users/tushevg/Desktop/dataTracks_17Apr2019.mat
%}
clearvars -except data metainfo

tShort = [15;30;45];
tLong = [90;120;150];
nFiles = size(metainfo, 1);
nGenes = size(data, 1);
k = 500;
thresh = 0.85;

idxUse = (metainfo(:,1) < 90) & (metainfo(:,2) < 4);

[~,~,idxMerge] = unique(metainfo(idxUse,1));


res = zeros(nGenes, 3);
plotData = zeros(4, 400, nGenes);
plotNorm = 0;

for k = 1 : nGenes
    
    nCodons = size(data(k).tracks, 2) / 3;
    
    if nCodons >= 400
        
        tracks = double(data(k).tracks(idxUse,:) + 1);
        fctr = mean(tracks(:,300:end-20), 2);
        fctr = mean(fctr) ./ fctr;
        tracks = bsxfun(@times, tracks, fctr);
        
        codons = zeros(4, nCodons);
        tShort = [0,15,30,45];
        tArray = [];
        xArray = [];
        for c = 1 : 4
            yN = sum(tracks(idxMerge == c,:), 1);
            yC = sum(reshape(yN, 3, nCodons, 1), 1);
            yC = movmean(yC, 30);
            %yC = yC ./ mean(yC(300:end-20));
            codons(c,:) = yC;
            
            if c > 1
                idxThresh = codons(c,:) > codons(1,:) .* thresh;
                slPoints = find(idxThresh(2:end) == 1 & idxThresh(1:end-1) == 0) + 1;
                slUse = false(length(slPoints),1);
                if ~isempty(slPoints)
                    for s = 1 : length(slPoints)
                        i = slPoints(s);
                        j = min([slPoints(s)+19,nCodons]);
                        slUse(s) = sum(idxThresh(i : j)) >= 10;
                    end
                else
                    slPoints(1) = nCodons;
                    slUse(1) = true;
                end
                
                tArray = [tArray;tShort(c)];
                xArray = [xArray;min(slPoints(slUse))];
                
            end
            
        end
        
        
        if issorted(xArray) && all(tArray == tShort(2:end)') && length(xArray) == 3
            p = polyfit(tArray, xArray, 1);
            res(k,1) = p(1);
            res(k,2) = nCodons;
            res(k,3) = corr(tArray, xArray).^2;
            if any(isnan(codons(1,1:400)))
                break;
            end
            plotData(1,:,k) = plotData(1,:,k) + codons(1,1:400);
            plotData(2,:,k) = plotData(2,:,k) + codons(2,1:400);
            plotData(3,:,k) = plotData(3,:,k) + codons(3,1:400);
            plotData(4,:,k) = plotData(4,:,k) + codons(4,1:400);
            plotNorm = plotNorm + 1;
        end
        
        
    end
    
end


tmp = [{data.gene}',num2cell(res)];
idxFinal = res(:,1) > 0;
tmp = tmp(idxFinal,:)';
fw = fopen('/Users/tushevg/Desktop/elongationTest_02May2019.txt','w');
fprintf(fw, '%s\t%.4f\t%d\t%.4f\n',tmp{:});
fclose(fw);




Z = mean(plotData(1,:,:), 3);
R = mean(plotData(2,:,:), 3);
G = mean(plotData(3,:,:), 3);
B = mean(plotData(4,:,:), 3);

figure('color','w');
hold on;
plot(thresh.*Z,'k');
plot(R,'r');
plot(G,'g');
plot(B,'b');
hold off;

%{
tmp = zeros(nGenes,2);

check = zeros(4,400);

for k = 1 : nGenes
    
    
    tracksShort = [sum(data(k).tracks(1:3,:));...
                   sum(data(k).tracks(13:15,:));...
                   sum(data(k).tracks(16:18,:));...
                   sum(data(k).tracks(19:21,:))];
    tracksLong = [sum(data(k).tracks(4:6,:));...
                  sum(data(k).tracks(22:24,:));...
                  sum(data(k).tracks(7:9,:));...
                  sum(data(k).tracks(10:12,:))];
            
    
    nCodons = size(data(k).tracks, 2) / 3;
    tmp(k,2) = nCodons;
    
    
    
    if (nCodons >= 400)
        
        %{
        [tSL, xSL] = getSLpoints(tracksShort, 400, tShort, 0.75);
        if issorted(xSL) && length(xSL) > 1
            p = polyfit(tSL, xSL, 1);
            tmp(k) = p(1);
            
        end
        %}
        
        res = zeros(4, nCodons);
        fctr = zeros(4, 1);
        
        for t = 1 : 4
            yN = tracksShort(t,:);
            yC = sum(reshape(yN, 3, nCodons, 1));
            yC = movmean(yC, 30);
            %yC = yC ./ mean(yC(600:end-20));
            fctr(t) = sum(yC(300:end-20)); 
            res(t,:) = yC;
        end
        fctr = mean(fctr)./fctr;
        res = bsxfun(@times, res, fctr);
        
        if ~any(any(isnan(res),2)) && ~any(any(isinf(res),2))
        check(1,:) = check(1,:) + res(1,1:400);
        check(2,:) = check(2,:) + res(2,1:400);
        check(3,:) = check(3,:) + res(3,1:400);
        check(4,:) = check(4,:) + res(4,1:400);
        i = i + 1;
        end
        
        %break;
        
    end
    
    
    %{
    if (nCodons >= 790)
        [tSL, xSL] = getSLpoints(tracksShort, 400, tLong, 0.8);
        %[tSL, xSL] = getSLpoints(tracksLong, 750, tLong, 0.5);
        if issorted(xSL) && length(xSL) > 1
            p = polyfit(tSL, xSL, 1);
            tmp(k) = p(1);
            
        end
        
    end
    %}
    
    
end

check = check ./ i;

figure('color','w');
hold on;
plot(check(1,:).*1, 'k');
plot(check(2,:), 'r');
plot(check(3,:), 'g');
plot(check(4,:), 'b');
hold off;

%}


%
function [tSL,xSL] = getSLpoints(tracks, xNorm, tNorm, thresh)

    [nTracks, nCodons] = size(tracks);
    nCodons = nCodons/3;
    
    % smooth and normalize
    res = zeros(nTracks, nCodons);
    xSL = zeros(nTracks - 1, 1);
    slUse = false(nTracks - 1, 1);
    for t = 1 : nTracks
        yN = tracks(t,:);
        yC = sum(reshape(yN, 3, nCodons, 1));
        yC = movmean(yC, 30);
        yC = yC ./ mean(yC(xNorm:end-20));
        res(t,:) = yC;
        
        if (t > 1)
            idx = res(t,:) > (res(1,:) .* thresh);
            
            slPoints = find((idx(1:end-1) == 0) & (idx(2:end) == 1)) + 1;
            
            if ~isempty(slPoints)
            
                idxUse = false(length(slPoints),1);
                for s = 1 : length(slPoints)
                    idxUse(s) = sum(idx(slPoints(s):min([slPoints(s)+19,nCodons]))) >= 10;
                end

                if sum(idxUse) > 0
                    xSL(t-1,1) = min(slPoints(idxUse));
                    slUse(t-1) = true;
                else
                    xSL(t-1,1) = nCodons;
                end
            
            end
            
            
        end
        
    end
    
    xSL = xSL(slUse);
    tSL = tNorm(slUse);
    
    
    

    
end

