%% plotPSites
clc
clear variables
close all

fh = fopen('/Users/tushevg/Desktop/tablePSitesProjection_NPL_MONO_12Feb2019.txt','r');
txt = textscan(fh, '%s %s %s','delimiter','\t');
fclose(fh);
fname = txt{1};
val = txt{2};
cnt = txt{3};

Y = zeros(9, 300);
idxStart = (1:100);
idxCenter = (101:200);
idxEnd = (201:300);
for k = 1 : length(val)
    
    v = sscanf(val{k}, '%f,');
    c = sscanf(cnt{k}, '%f,');
    Y(k,idxStart) = v(idxStart)./mean(v(idxStart));
    Y(k,idxCenter) = v(idxCenter)./mean(v(idxCenter));
    Y(k,idxEnd) = v(idxEnd)./mean(v(idxEnd));
    
end

%
idx = 1:3;
titleSample = 'Neuropil-Poly';

idxX = [-25:74;...
        81:180;...
        186:285];
idxY = [1:100;...
        101:200;...
        201:300];
  
clr_mtx = [30,144,255;...
           148,0,211;...
           144,238,144]./255;

figure('color','w');
hold on;
for k = 1 : 3
    
    Ym = mean(Y(idx,idxY(k,:)), 1);
    Yse = std(Y(idx,idxY(k,:)), [], 1);
    Xbin = idxX(k,:);
    xpatch = [Xbin, fliplr(Xbin)];
    ypatch = [(Ym-Yse), fliplr((Ym+Yse))];
    
    hp = fill(xpatch, ypatch, clr_mtx(k,:));
    set(hp,'edgecolor','none','facealpha',0.5);
    plot(Xbin, Ym, 'color',clr_mtx(k,:));
    
end
hold off;

xtick = [-25,0,25,50,...
         81,106,131,156,180,...
         211,236,261,285];
xticklabel = num2cell([(-25:25:50),...
              (-50:25:50),...
              (-50:25:26)]);
xticklabel(2) = {'start'};
xticklabel(7) = {'center'};
xticklabel(12) = {'stop'};
          
set(gca,'xtick',xtick,...
        'xticklabel',xticklabel,...
        'ylim',[0,8]);
xlabel('relative offset [nts]','fontsize',12);
ylabel('relative coverage of p-site','fontsize',12);
title(titleSample,'fontsize',12,'fontweight','normal');
print(gcf, '-dsvg','-r300',sprintf('figurePSiteProjection_%s.svg',titleSample));

nfft = 256; % next larger power of 2
Fs = 1;

Xff = (0:nfft/2)* Fs/nfft;
Yff = zeros(length(idx),ceil(nfft/2+1));
for k = 1 : length(idx)
    Yn = mean(Y(idx(k),:), 1);
    Yz = Yn - mean(Yn);

    Yf = fft(Yz,nfft); % Fast Fourier Transform
    Yp = abs(Yf.^2); % raw power spectrum density
    
    Yff(k,:) = Yp(1:1+nfft/2); % half-spectrum
end

Ym = mean(Yff,1);
Yse = std(Yff,[],1)./sqrt(size(Yff,1));
xpatch = [Xff, fliplr(Xff)];
ypatch = [(Ym-Yse), fliplr((Ym+Yse))];
Xexp = 1/3;
figure('color','w');
hold on;
hp = fill(xpatch, ypatch, [0.5,0.5,0.5]);
set(hp,'edgecolor','none','facealpha',0.5);
plot(Xff, Ym, 'color','k');
ylim = get(gca,'ylim');
plot([Xexp,Xexp],[0,0.025*max(ylim)],'r','linewidth',1.2);
hold off;
ylabel('power spectrum desnity','fontsize',12);
xlabel('frequency [Hz]');
title(titleSample,'fontsize',12,'fontweight','normal');
print(gcf, '-dsvg','-r300',sprintf('figurePSiteFTransform_%s.svg',titleSample));