clear all;
close all;
clc;

sourcepath = {'C:\Users\Delmon\Dropbox\Jobb\Matlab\Signal denoising\AIPh.xlsx', ...
    'C:\Users\Delmon\Dropbox\Jobb\Matlab\Signal denoising\Amphinex.xlsx', ...
    'C:\Users\Delmon\Dropbox\Jobb\Matlab\Signal denoising\Temoporfin.xlsx', ...
    'C:\Users\Delmon\Dropbox\Jobb\Matlab\Signal denoising\spectrum1.xlsx', ...
    'C:\Users\Delmon\Dropbox\Jobb\Matlab\Signal denoising\spectrum2.xlsx'};

minvec = [0 1000 1000 1000 1000];
maxvec = [950 1750 1950 1750 1750];
SIvalues = [];

for i = 4 % 1:length(sourcepath)
    [SI, lambda] = spectrumSmoothing(sourcepath{i}, minvec(i), maxvec(i));
%     [~, ind] = min(abs(lambda - 700)); % find(lambda == 700);
%     SIvalues = [SIvalues; SI(ind)];
%     struct{i}.SI = SI;
%     struct{i}.lambda = lambda;
end

% concentations = [10^(-3) 10^(-4) 10^(-5)];
% [r, pval] = corr(SIvalues, concentations, 'type', 'Pearson');
% X = [ones(length(concentations),1) concentations];
% beta = X\log(SI);

