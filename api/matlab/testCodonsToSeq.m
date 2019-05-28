%% testCodonsToSeq
clc
close all

%{

clear variables
fh = fopen('testCoverage.txt', 'r');
txt = textscan(fh, '%s %n', 'delimiter', '\t');
fclose(fh);
S = txt{1};
D = txt{2};

S(isinf(D)) = [];
D(isinf(D)) = [];
[codonList, ~, codonIndex] = unique(S);
%}
clearvars -except S D codonList codonIndex;

% x(1) = m
% x(2) = sigma
% x(3) = lambda
yf = @(b, x) (b(3)/2).*exp((b(3)./2).*(2.*b(1)+b(3).*b(2)^2 - 2.*x)).*erfc((b(1)+b(3).*b(2).^2 - x) ./ (sqrt(2).*b(2)));

xbin = linspace(0,10,100)';

codonFreq = zeros(64, 1);
tic
k = 1;
for k = 1 : 64
    Y = D(codonIndex == k);
    Y(Y>10) = [];
    
    lambda = std(Y) .* 0.8;
    mu = mean(Y) - skewness(Y);
    sigma = sqrt(var(Y) - (lambda^2));
    b0 = [mu, sigma, lambda];
    
    [Yfit, Xfit] = hist(Y, 128);
    
    mdl = fminsearch(@(params)logLHexGauss(params, Y), b0);
    codonFreq(k) = mdl(1);
    
    
    %{
    Ytmp = exGaussPDF(mdl, Xfit);
    figure('color','w')
    hold on;
    plot(Xfit, Yfit./max(Yfit), 'k.','markersize',10);
    plot(Xfit, Ytmp./max(Ytmp), 'r','linewidth',1.2);
    plot([mdl(1),mdl(1)],[0,1],'g','linewidth',1.2);
    
    hold off;
    %}
    
end
toc

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

