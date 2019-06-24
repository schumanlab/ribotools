%% calculateNFCperList
clc
clear variables
close all

%% readCodonTable
table = readPlainCodonTable('/Users/tushevg/Desktop/data/codonTable.txt');

%% read codons list
fh = fopen('codonsNFC_NplTotal_18Jun2019.txt','r');
txt = textscan(fh,'%s %n','delimiter','\t');
fclose(fh);
codon = txt{1};
value = txt{2};

%% calculate parameters
a = 5;

params = zeros(64,3);

for a = 1 : 64
    
    idx_aa = strcmp(table.codon{a}, codon);
    Y = value(idx_aa);
    lambda = std(Y) .* 0.8;
    mu = mean(Y) - skewness(Y);
    sigma = sqrt(var(Y) - (lambda^2));
    b0 = [mu, sigma, lambda];

    params(a,:) = fminsearch(@(params)logLHexGauss(params, Y), b0);
    
end


ATR = 4;
timeDecoding = params(:,1);
%timePausing = params(:,1) + 3.* params(:,2) + abs(log(0.5))./params(:,3);
timePausing = abs(log(0.5))./params(:,3);

timeDecoding = (1/ATR) .* (timeDecoding./mean(timeDecoding));
timePausing = (1/ATR) .* (timePausing./mean(timeDecoding));
X = exGaussPDF(params, timeDecoding);


table.w = [timeDecoding, timePausing];
writeCodonTable(table, '/Users/tushevg/Desktop/codonsParams_NplTotal_18Jun2019.txt');

function table = readPlainCodonTable(fileName)
    
    fh = fopen(fileName, 'r');
    txt = textscan(fh,'%s %s %s %s','delimiter','\t');
    fclose(fh);
    table.codon = txt{1};
    table.letter = txt{2};
    table.code = txt{3};
    table.name = txt{4};
end



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

function writeCodonTable(table, fileName)

    out = [table.codon, table.letter, table.code, table.name, num2cell(table.w)]';
    
    fw = fopen(fileName,'w');
    fprintf(fw, '%s\t%c\t%s\t%s\t0\t%.8f\t%.8f\n',out{:});
    fclose(fw);

end