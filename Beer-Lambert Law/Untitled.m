clear all;
close all
clc

path = 'C:\Users\Delmon\Dropbox\Jobb\Matlab\Beer-Lamberts Law';
sourcefile = 'PpIX_Ppp.xlsx';

[c, cint, stats] = EstimateConcentrationsBeerLambertLaw( ... 
    fullfile(path, sourcefile));

