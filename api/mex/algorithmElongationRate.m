%% algorithmElongationRate
clc
close all

%% load tracks
%{
clear variables
load /Users/tushevg/Desktop/data/timepoints/dataTracks_17Apr2019.mat
%}
clearvars -except data metainfo

timeList = metainfo(:,1);
timeList(4:6) = 1;
[timeListUnique,~,idxMerge] = unique(timeList);


nTimes = max(idxMerge);
nGenes = length(data);
g = 2; % 8359 longest


arrayTime = [15;30;45;90;120;150];

%for q = 3 : 5
nBin = q;
X = arrayTime(1:nBin);
Xfit = [ones(nBin,1), X];
result = zeros(nGenes, 3);
tic
%for g = 1 : nGenes

    tracks = double(data(g).tracks);
    %counts = double(data(g).counts);
    %factor = geomean(counts)./counts;
    %tracks = bsxfun(@times, tracks, factor);
    nCodons = size(tracks, 2) / 3;
    codons = zeros(nTimes, nCodons);
    for t = 1 : nTimes
        arrayNucl = sum(tracks(idxMerge == t,:) + 1, 1);
        arrayCodons = sum(reshape(arrayNucl, 3, nCodons), 1);
        arrayCodons = movmean(arrayCodons, 5);
        arrayCodons = arrayCodons ./ mean(arrayCodons(end-ceil(0.1*nCodons):end));
        codons(t,:) = arrayCodons;
    end
    
    %% signal squered error
    sse = [(codons(3,:) - codons(1,:));...
           (codons(4,:) - codons(1,:));...
           (codons(5,:) - codons(1,:));...
           (codons(6,:) - codons(2,:));...
    	   (codons(7,:) - codons(2,:));...
           (codons(8,:) - codons(2,:))];
    %}
    %sse = bsxfun(@rdivide, sse, std(codons(3:8,:), [] ,2));
    
    %% threshold crossing
    Xraw = [];
    Yraw = [];
    for j = 1 : nBin
        sseNow = sse(j,:);
        crossNow = find((sseNow(1:end - 1) < 0) & (sseNow(2:end) >= 0));
        crossNow(crossNow > X(j) * 10) = [];
        Xraw = cat(1, Xraw, repmat(arrayTime(j),length(crossNow), 1));
        Yraw = cat(1, Yraw, crossNow');
    end
    
    [~, ~, idxGroup] = unique(Xraw);
    simSize = histc(idxGroup, (1:nBin));
    n = max(idxGroup);
    if isempty(Xraw)
        n = 0;
    end
    
    if n == nBin
        
       if nBin == 5
            D = combvec(Yraw(idxGroup == 1)',...
                        Yraw(idxGroup == 2)',...
                        Yraw(idxGroup == 3)',...
                        Yraw(idxGroup == 4)',...
                        Yraw(idxGroup == 5)');
       end

       if nBin == 4
            D = combvec(Yraw(idxGroup == 1)',...
                        Yraw(idxGroup == 2)',...
                        Yraw(idxGroup == 3)',...
                        Yraw(idxGroup == 4)');
       end

       if nBin == 3
            D = combvec(Yraw(idxGroup == 1)',...
                        Yraw(idxGroup == 2)',...
                        Yraw(idxGroup == 3)');
       end
       
       rsq = zeros(size(D,2),1);
       for j = 1 : size(D,2)
           y = D(:,j);
           b = Xfit \ y;
           
           if (b(1) > 0)
            c = corrcoef(X, y);
            rsq(j,1) = 1 - 4/3 * (1 - c(1,2));
           end
           
       end
       [rsqmax,idx] = max(rsq);
       p = robustfit(X,D(:,idx));
        
       result(g,1) = nCodons;
       result(g,2) = p(2);
       result(g,3) = rsqmax;
        
        
    end
    %
    
    %Dx = pdist(Xraw, @(x,y)(y-x));
    %Dy = pdist(Yraw, @(x,y)(y-x));
    %tmp = ones(size(Xraw,1));
    %tmp = tril(tmp, -1);
    %[rowIdx, colIdx] = find(tmp);
    
    %Dy(Dx==0) = [];
    %Dx(Dx==0) = [];
    %D = Dy ./ Dx;
    
    %{
    if any(simSize == 0)
        %continue;
    end
    
    [i,j,l,m,n] = ndgrid(Yraw(idxGroup==1),...
                         Yraw(idxGroup==2),...
                         Yraw(idxGroup==3),...
                         Yraw(idxGroup==4),...
                         Yraw(idxGroup==5));
    D = [i(:),j(:),l(:),m(:),n(:)];
    dD = bsxfun(@rdivide, diff(D,[],2), dX);
    iD = sum(dD, 2);
    [~,idx] = min(iD);
    
    
    p = robustfit(X,D(idx,:)');
    
    
    result(g,1) = nCodons;
    result(g,2) = p(2);
    result(g,3) = corr(X,D(idx,:)').^2;
    Yfit = polyval(flipud(p), X);
    %}
    %{
    figure('color','w');
    plot(Xraw, Yraw, '.', 'color', [.65,.65,.65],'markersize',10);
    hold on;
    plot(X,D(:,idx),'k.','markersize',15);
    %plot(X, Yfit,'k');
    hold off;
    set(gca,'box','off',...
            'xlim',[0,120],...
            'xtick',X,...
            'ylim',[0,max(Yraw)]);
    %}
    
    
    
%end
%fprintf('Finished %d in %.8f\n',nBin, toc);
%fileName = sprintf('/Users/tushevg/Desktop/elongationRate_TP%d_16May2019.mat',nBin);
%save(fileName, 'result');

%end

