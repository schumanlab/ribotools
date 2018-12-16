%% plotPSites
clc
clear variables
close all

fh = fopen('tablePsiteProject_MonoVsPoly.txt','r');
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
    c = sscanf(cnt{k}, '%d,');
    Y(k,idxStart) = v(idxStart)./mean(c(idxStart));
    Y(k,idxCenter) = v(idxCenter)./mean(c(idxCenter));
    Y(k,idxEnd) = v(idxEnd)./mean(c(idxEnd));
    
end

idx = 7:9;
Xbin = (1:300)';
Ym = mean(Y(idx,:),1);
Yse = std(Y(idx,:),[],1);

figure('color','w');
hold on;
xbin = -25:74;
plot(xbin,Ym(1:100));

xbin = 81:180;
plot(xbin, Ym(101:200));

xbin = 186:285;
plot(xbin, Ym(201:300));

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
        'xticklabel',xticklabel);


%{
xpatch = [Xbin', fliplr(Xbin')];
ypatch = [(Ym-Yse), fliplr((Ym+Yse))];

figure('color','w');
hold(gca,'on');
hp = fill(xpatch, ypatch, [.75, .75, .75]);
plot(Xbin, Ym, 'k', 'linewidth', 1);
set(hp, 'EdgeColor', 'none');
hold(gca,'off');
%}

%{
Yz = Yn - mean(Yn);

nfft = 256; % next larger power of 2
Yf = fft(Yz,nfft); % Fast Fourier Transform
Yp = abs(Yf.^2); % raw power spectrum density
Fs = 1;
f_scale = (0:nfft/2)* Fs/nfft;
Yhs = Yp(1:1+nfft/2); % half-spectrum

figure('color','w');
plot(Yn,'k','linewidth',1.2);

figure('color','w');
plot(f_scale, Yhs, 'k','linewidth',1.2);
%}