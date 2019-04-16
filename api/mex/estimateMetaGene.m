%% estimateMetaGene
clc
close all

%{
clear variables
load dataTracks_TimePoints.mat
%}
clearvars -except bed metainfo;

k = 703;
hFilter = ones(1, 5);
nfiles = length(metainfo);
idx = zeros(length(bed),1);

idxshort = metainfo(:,2) <= 3;
idxlong = metainfo(:,2) > 3;
bin = 800;
data = zeros(1, bin, nfiles);
i = 0;

for k = 1 : length(bed)

    cdsStart = bed(k).cdsStart;
    cdsEnd = bed(k).cdsEnd;
    cdsSpan = cdsEnd - cdsStart;
    offset = mod(cdsSpan, 3);
    
    cdsEnd = cdsEnd - offset;
    cdsSpan = cdsEnd - cdsStart;
    %
    
    track = bed(k).linearCoverage;
    track = double(track(:,cdsStart + 1 : cdsEnd));
    
    %{
    i = i + 1;
    c = zeros(nfiles,1);
    f = 1;
    for f = 1 : nfiles
        trackCodons = track(f,:);
        trackCodons = sum(reshape(trackCodons, 3, length(trackCodons)/3), 1);
        trackCodons = conv(trackCodons, hFilter, 'same');
        trackCodons = trackCodons(1:bin) + 1;
        
        %x = (1:length(trackCodons));
        trackCodons = trackCodons ./ mean(trackCodons(end-200:end));
        %c(f, 1) = find((trackCodons > 0.5) & (x > 40), 1, 'first');
        data(i, :, f) = trackCodons;
    end
    %}
  
    
    %{
    figure('color','w');
    hold on;
    plot(metainfo(idxshort,1), c(idxshort), 'k.');
    plot(metainfo(idxlong,1), c(idxlong), 'r.');
    
    hold off;
    %}
    
    
    
    %
    if (cdsSpan < bin*3)
        continue;
    end
    
    i = i + 1;
    
    
    track = bed(k).linearCoverage;
    track = track(:,cdsStart + 1 : cdsEnd);
    for f = 1 : nfiles
        trackCodons = track(f,:);
        trackCodons = sum(reshape(trackCodons, 3, length(trackCodons)/3), 1);
        trackCodons = conv(trackCodons, hFilter, 'same');
        trackCodons = trackCodons(1:bin) + 1;
        trackCodons = trackCodons ./ mean(trackCodons(end-100:end));
        data(i,:, f) = trackCodons;
    end
    %}
end

%% normalize to 0-time point
baseShort = squeeze(median(data(:,:,1:3), 3));
baseLong = squeeze(median(data(:,:,4:6), 3));
idxshort = metainfo(:,2) <= 3;
idxlong = metainfo(:,2) > 3;
dataShort = bsxfun(@rdivide, data(:,:,idxshort), baseShort);
dataLong = bsxfun(@rdivide, data(:,:,idxlong), baseLong);

data0 = median(mean(cat(3, dataShort(:,:,1:3), dataLong(:,:,1:3)), 3), 1);
data15 = median(mean(dataShort(:,:,4:6), 3), 1);
data30 = median(mean(dataShort(:,:,7:9), 3), 1);
data45 = median(mean(dataShort(:,:,10:12), 3), 1);
data90 = median(mean(dataLong(:,:,4:6), 3), 1);
data120 = median(mean(dataLong(:,:,7:9), 3), 1);
data150 = median(mean(dataLong(:,:,10:12), 3), 1);

thresh_x = 50;
thresh_y = 0.85;%data15(thresh_x);

figure('color','w');
hold on;
plot([0,800],[1,1],'color',[.25,.25,.25]);
plot([0,800],[thresh_y, thresh_y],'color',[.25,.25,.25]);
plot([thresh_x,thresh_x],[0,1.5], 'color', [.25, .25, .25]);
h(1) = plot(data0,'linewidth',1.2);
h(2) = plot(data15,'linewidth',1.2);
h(3) = plot(data30,'linewidth',1.2);
h(4) = plot(data45,'linewidth',1.2);
h(5) = plot(data90,'linewidth',1.2);
h(6) = plot(data120,'linewidth',1.2);
h(7) = plot(data150,'linewidth',1.2);

hold off;
hl = legend(h, '0 sec', '15 sec', '30 sec', '45 sec', '90 sec', '120 sec', '150 sec');
set(hl, 'edgecolor', 'w', 'location', 'southeast', 'fontsize', 14);
set(gca,'fontsize',14);
xlabel('codons','fontsize',14);
ylabel('normalized coverage','fontsize',14);
print(gcf, '-dsvg', '-r300', 'fiture_AvgTranslationRateMetaGene_12Apr2019.svg');


x = (1:800);
y = zeros(6,1);
y(1) = find(data15 >= thresh_y & x >= thresh_x, 1, 'first');
y(2) = find(data30 >= thresh_y & x >= thresh_x, 1, 'first');
y(3) = find(data45 >= thresh_y & x >= thresh_x, 1, 'first');
y(4) = find(data90 >= thresh_y & x >= thresh_x, 1, 'first');
y(5) = find(data120 >= thresh_y & x >= thresh_x, 1, 'first');
y(6) = find(data150 >= thresh_y & x >= thresh_x, 1, 'first');

t = unique(metainfo(:,1));
t(1) = [];
p = polyfit(t, y, 1);
tfit = linspace(5,t(end), 100);
yfit = polyval(p, tfit);
figure('color','w');
hold on;
plot(tfit, yfit, 'color', [.65,.65,.65]);
plot(t, y, 'k.', 'markersize', 12);
hold off;
set(gca,'box','off','xlim',[0,150],'xtick',t,'ylim',[0,450],'fontsize',14);
xlabel('time [sec]','fontsize',14);
ylabel('crossing point [codons]','fontsize',14);
text(15, 400, sprintf('y(t) = %.4f * t - %.4f\nR^2 = %.4f',p(1), abs(p(2)), corr(t,y).^2),'fontsize',14);
print(gcf, '-dsvg', '-r300', 'fiture_AvgTranslationRateSlope_12Apr2019.svg');


%{
dt = metainfo(:,1);
dt(4:6) = 1;
tmp = squeeze(mean(data, 1));
tu = unique(dt);
mtx = zeros(bin, length(tu));
for t = 1 : length(tu)
    idx = tu(t) == dt;
    mtx(:,t) = mean(tmp(:,idx),2);
    
end

%mtx(:,1) = mtx(:,1) * 0.8;
%mtx(:,2) = mtx(:,2) * 0.5;

A = bsxfun(@rdivide, mtx(:,[1,3,4,5]), mtx(:,1));
B = bsxfun(@rdivide, mtx(:,[2,6,7,8]), mtx(:,2));



figure('color','w');
plot([A,B]);
%}

%% calculate average slope
%{
thresh = linspace(0.4, 0.9, 100);
crr = zeros(length(thresh), 1);
slp = zeros(length(thresh), 1);
for j = 1 : length(thresh)

c = zeros(length(tu),1);
x = (1:bin)';
for k = 1 : size(mtx, 2)
    
    c(k, 1) = find((mtx(:,k) > thresh(j)) & (x > 15), 1, 'first');
    
end

crr(j) = corr(tu, c);
p = polyfit(tu, c, 1);
slp(j) = p(1);



end

%% figure
figure('color','w');
plot(thresh, slp);

figure('color','w');
plot(thresh, crr);

%% figure
idxshort = [1,3,4,5];
idxlong = [2,6,7,8];
figure('color','w')
plot(tu, c, 'k.');

pshort = polyfit(tu(idxshort), c(idxshort), 1);
plong = polyfit(tu(idxlong), c(idxlong), 1);
p = polyfit(tu, c, 1);

%}

%% read bed file
%{
bed = readBedFile('/Users/tushevg/Desktop/data/bed/ncbiRefSeq_rn6_IP_29Jan2019_appris.bed');
name = {bed.transcript}';
cdsStart = [bed.cdsStart]';
cdsEnd = [bed.cdsEnd]';
cdsSpan = [bed.cdsSpan]';
offset = mod(cdsSpan, 3);
idxOne = (offset == 1);
idxTwo = (offset == 2);
offset(idxOne) = 2;
offset(idxTwo) = 1;
cdsEnd = cdsEnd + offset;

%% fileList
fileList = ls('/Users/tushevg/Desktop/data/timepoints/*.gbed.gz');
fileList(end) = [];
fileList = regexp(fileList, '\n', 'split')';
nFiles = length(fileList);
fileName = regexprep(fileList, '/Users/tushevg/Desktop/data/timepoints/','');
fileName = regexprep(fileName, '_transcriptome_sorted_umi.gbed.gz', ''); 
fileName = cellfun(@(x) {sscanf(x,'%ds_%d')'},fileName);
metainfo = cell2mat(fileName);

%% allocate handels
handles = zeros(nFiles, 1, 'uint64');
for k = 1 : nFiles
    handles(k) = coverage('new', fileList{k});
end

%% track matrix
txName = 'NM_012920.1';
idx = strcmp(txName, name);
k = 1;
tic
for k = 1 : length(bed)
    
    tracks = zeros(nFiles, bed(k).geneSpan, 'int32');
    
    for f = 1 : nFiles
        tracks(f,:) = coverage('query', handles(f), bed(k).transcript, 0, bed(k).geneSpan);
    end
    
    bed(k).linearCoverage = tracks;

end
toc
save('dataTracks_TimePoints.mat','bed', 'metainfo')
%}


%{
txStart = 0;
txEnd = 1498;
txSpan = txEnd - txStart;

%}


%{
%% deallocate handels
for k = 1 : nFiles
    coverage('delete', handles(k));
end

%}






%% --- FUNCTIONS --- %%
function bed = readBedFile(fileBed)

    fh = fopen(fileBed, 'r');
    txt = textscan(fh, '%s', 'delimiter', '\n');
    fclose(fh);
    
    raw = txt{1};
    bed = repmat(BedLine(), length(raw), 1);
    
    for b = 1 : length(raw)
        tmp = BedLine();
        tmp.parse(raw{b});
        bed(b) = tmp;
    end

end