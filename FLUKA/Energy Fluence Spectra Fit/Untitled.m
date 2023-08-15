clear all;
close all;
fclose('all');
clc;

path = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\X-ray 220 keV\Energy Fluence Spectra';
sourcefilename = {'USRBDX1.xlsx'};
destfilename = {'Xray220kVFluenceEnergyDoubleGaussFit.xlsx'};
 
% str = {{'R_{40}', 'R_{60}', 'R_{80}', 'BP MC', ...
%     'Distal R_{80}', 'Distal R_{60}', 'Distal R_{40}'}};
str = {{'220 kV X-ray, 1.52 mm Al + 0.70 mm Cu filter', ...
    '220 kV X-ray, no filter'}};

ndetectors = [2];
nbins = [400];

for i = 1
    EnergyFluenceGaussianFit(path, sourcefilename{i}, destfilename{i}, ...
        ndetectors(i), nbins(i), str{i})
end