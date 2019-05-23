%% clusteringMonoVsPoly
clc
close all

%{
clear variables

%% read ibaq
fh = fopen('/Users/tushevg/Desktop/data/clustering/20190223_compartments_ibaq_matrix.txt', 'r');
hd = fgetl(fh);
txt = textscan(fh,'%n %n %n %n %n %n %s', 'delimiter', '\t');
fclose(fh);
ibaq.name = txt{7};
ibaq.counts = [txt{1:end-1}];
ibaq.npl = nanmean(ibaq.counts(:,1:3), 2);
ibaq.smt = nanmean(ibaq.counts(:,4:6), 2);

idxFilter = cellfun('isempty',ibaq.name) | isnan(ibaq.npl) | isnan(ibaq.smt);
ibaq.name(idxFilter) = [];
ibaq.npl(idxFilter,:) = [];
ibaq.smt(idxFilter,:) = [];
[~,idxUnq] = unique(ibaq.name);

ibaq.name = ibaq.name(idxUnq);
ibaq.npl = ibaq.npl(idxUnq);
ibaq.smt = ibaq.smt(idxUnq);
nibaq = length(ibaq.name);



%% read mono/poly footprint
fh = fopen('/Users/tushevg/Desktop/data/clustering/20190302_neuropil_02_df_elongating.txt', 'r');
hd = fgetl(fh);
txt = textscan(fh, '%s %n %n %n %s %s', 'delimiter', '\t');
fclose(fh);
fpr.name = txt{1};
fpr.exp = txt{2};
fpr.log2FC = txt{3};
fpr.padj = txt{4};
fpr.enrichment = txt{5};
fpr.celltype = txt{6};
idxCellType = strcmp('neuronal',fpr.celltype) & (fpr.padj > 0) & (fpr.padj <= 0.05);


%% read rna
fh = fopen('/Users/tushevg/Desktop/data/clustering/totalRNA_exon_tpms_07May.txt', 'r');
hd = fgetl(fh);
txt = textscan(fh, '%s %*n  %*n  %*n  %*n %n', 'delimiter', '\t');
fclose(fh);
rna.name = txt{1};
rna.npl = txt{2};

%% merage master table
table = zeros(nibaq, 4);
idxfilter = false(nibaq,1);
for k = 1 : nibaq
    idxMP = strcmp(ibaq.name{k}, fpr.name) & idxCellType;
    idxRN = strcmp(ibaq.name{k}, rna.name);
    if sum(idxMP) > 0
        table(k,:) = [fpr.exp(idxMP), fpr.log2FC(idxMP), rna.npl(idxRN), ibaq.npl(k)];
        idxfilter(k) = true;
    end
end
table(~idxfilter,:) = [];
table(:,4) = table(:,4) - 8.4;
names = ibaq.name(idxfilter);
%}
clearvars -except table names;

%% clustering
mtx = zscore([table(:,2), log2(table(:,3)), table(:,4)]);

D = pdist(mtx, 'euclidean');
Z = linkage(D, 'complete');
idx = cluster(Z, 'maxclust', 16);


out = [names, num2cell(table), num2cell(mtx)]';
fh = fopen('tableX_RNAvsProtein_13May2019.txt','w');
fprintf(fh, 'name\tmean.expression\tfold.change\tRNA[TPM]\tProtein[IBAQ]\tzscore.FC\tzscore.RNA\tzscore.Protein\n');
fprintf(fh, '%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',out{:});
fclose(fh);

%{
idx = mtx(:,1) < 0;

cMono = [30,144,255]./255;
cPoly = [255,69,0]./255;

figure('color','w');
hold on;
plot([-3,3], [0,0], 'k');
plot([0,0], [-3,3], 'k');
h(1) = plot(mtx(idx,2),mtx(idx,3),'.','color',cPoly);
h(2) = plot(mtx(~idx,2),mtx(~idx,3),'.','color',cMono);

pmono = polyfit(mtx(~idx,2),mtx(~idx,3),1);
ppoly = polyfit(mtx(idx,2),mtx(idx,3),1);

xfit = linspace(-3,3,100)';
ymono = polyval(pmono, xfit);
ypoly = polyval(ppoly, xfit);
plot(xfit, ymono,'color',cMono,'linewidth',1.2);
plot(xfit, ypoly,'color',cPoly,'linewidth',1.2);
hold off;
ylabel('protein [zscore]','fontsize',14);
xlabel('rna [zscore]','fontsize',14);
set(gca,'xlim',[-3.5,3.5]);
set(gca,'ylim',[-3.5,3.5]);
set(gca,'fontsize',12);
hl = legend(h,'poly','mono');
set(hl,'edgecolor','w','location','northwest','fontsize',14);
print(gcf,'-dsvg','-r300','figure_MonoVsPoly_RNAvsProtein.svg');

%}




%{
cgobj = clustergram(mtx,...
                    'Standardize', 'none',...
                    'Cluster', 'column',...
                    'RowPDist','euclidean',...
                    'Linkage', 'complete',...
                    'OptimalLeafOrder', true,...
                    'DisplayRange', 2,...
                    'ColumnLabels',{'Mono/Poly','RNA','Protein'},...
                    'ColumnLabelsRotate',0,...
                    'Dendrogram',2.5);
%}

%idx = kmeans(mtx, 5,'Distance','sqeuclidean');
%[~,idxsort] = sortrows([idx,mean(mtx,2)],[1,2]);
%figure('color','w');
%imagesc(mtx(idxsort,:),[-3,3]);


%{
[~,f] = pca(mtx);
figure('color','w');
plot(f(:,2),f(:,3),'k.');

%
T = tsne(mtx);
figure('color','w');
plot(T(:,1), T(:,2),'k.');
%}









