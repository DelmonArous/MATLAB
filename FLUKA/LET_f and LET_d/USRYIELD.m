clear all;
close all;
fclose('all');
clc;

%%

destpath = 'C:\Users\delmo\Desktop\AM paper';

x1  = 81.287;
x2  = 81.4;
n   = 70; 

pos_vec     = linspace(x1, x2, n);
iter_vec    = (1:length(pos_vec)-1);

unitBIN_start = 25;

%% Write to file

% XYP LETpla1    81.2886
% XYP LETpla2    81.2903

% TARGET1      5 +target +LETpla1
% TARGET2      5 +target -LETpla1 +LETpla2
% TARGET3      5 +target -LETpla2

% USRYIELD         123  ALL-PART      -26.   TARGET1   TARGET2          LET2
% USRYIELD        100.       0.1     1000. $AbsEnrgy        0.      2703 &

fid = fopen(fullfile(destpath, 'USRYIELD.txt'), 'wt');
fprintf(fid, 'XYP LETpla%i    %.4f\n', [iter_vec; pos_vec(2:end)]);
fprintf(fid, 'TARGET%i      5 +target -LETpla%i +LETpla%i\n', ...
    [iter_vec+1; iter_vec; iter_vec+1]);
fprintf(fid, ['USRYIELD         123  ALL-PART      -%i.   TARGET%i   TARGET%i          LET%i\n', ...
    'USRYIELD        100.       0.1     1000. $AbsEnrgy        0.      2703 &\n'], ...
    [iter_vec+unitBIN_start; iter_vec; iter_vec+1; iter_vec+1]);
fclose(fid);


