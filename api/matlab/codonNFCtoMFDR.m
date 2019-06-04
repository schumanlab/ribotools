%% codonNFCtoMFDR
clc
clear variables
close all

%% read reference codons
fh = fopen('/Users/tushevg/Desktop/ribotools/api/cpp/codonTable.txt', 'r');
txt = textscan(fh, '%s %s %s %s','delimiter','\t');
fclose(fh);
ref.code = txt{1};
ref.letter = txt{2};
ref.abbr = txt{3};
ref.name = txt{4};


%% parse files
pathToNFC = '/Users/tushevg/Desktop/nfc/';
fileList = dir([pathToNFC, filesep, '*.txt']);
fileCount = length(fileList);
f = 1;

mfdr = zeros(length(ref.code), fileCount);

fileTagList = repmat({'unknown'}, fileCount, 1);

for f = 1 : fileCount
    fileName = [pathToNFC, filesep, fileList(f).name];
    [~, fileTag] = fileparts(fileName);
    fileTag = regexprep(fileTag, 'codonNFC_','');
    fileTagList(f) = {fileTag};
    
    
    fh = fopen(fileName, 'r');
    txt = textscan(fh, '%s %n','delimiter','\t');
    fclose(fh);
    qry.code = txt{1};
    qry.value = txt{2};
    
    idxFilter = isinf(qry.value) | isnan(qry.value);
    qry.code(idxFilter) = [];
    qry.value(idxFilter) = [];
    
    for a = 1 : length(ref.code)
        idx = strcmp(ref.code{a}, qry.code);
        Y = qry.value(idx);
    
        lambda = std(Y) .* 0.8;
        mu = mean(Y) - skewness(Y);
        sigma = sqrt(var(Y) - (lambda^2));
        b0 = [mu, sigma, lambda];
    
        mdl = fminsearch(@(params)logLHexGauss(params, Y), b0);
        mfdr(a, f) = mdl(1);
        
    end
    %}
end

fw = fopen('/Users/tushevg/Desktop/codonNFC_Report_NplSmt_04Jun2019.txt','w');
fprintf(fw,'#code\tletter\tabbreviation\tname\t');
header = sprintf('%s\t',fileTagList{:});
header(end) = [];
fprintf(fw,'%s\n',header);
for k = 1 : length(ref.code)
    fprintf(fw,'%s\t',ref.code{k});
    fprintf(fw,'%s\t',ref.letter{k});
    fprintf(fw,'%s\t',ref.abbr{k});
    fprintf(fw,'%s\t',ref.name{k});
    data = sprintf('%.6f\t',mfdr(k,:));
    data(end) = [];
    fprintf(fw,'%s\n',data);
end
fclose(fw);


function logL = logLHexGauss(params, x)

    y = exGaussPDF(params, x) + eps;
    logL = -sum(log(y));

end

function f = exGaussPDF(params, x)
    
    mu = params(1);
    sigma = params(2);
    lambda = params(3);
    
    argGauss = exp((lambda ./ 2) .* (2 .* mu + lambda .* sigma .^ 2 - 2 .* x));
    argExp = erfc((mu + lambda .* sigma .^ 2 - x) ./ (sqrt(2) .* sigma));
    f = (lambda ./ 2) .* argGauss .* argExp;
    
end
