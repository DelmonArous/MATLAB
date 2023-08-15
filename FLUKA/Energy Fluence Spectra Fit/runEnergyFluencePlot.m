clear all;
close all;
fclose('all');
clc;

sourcepath = 'C:\Users\delmo\Desktop\AM paper\Energy Fluence Spectra\Proton 15.2 MeV Water target\XLSX';
filename = 'AMproton_v3_28_tab_81point3cmDist.xlsx';

% str = {'1.52 mm Al + 0.70 mm Cu filter', 'No filter'};
str = {'Entrance', 'Position 1', 'Position 2', 'Position 3', ...
    'Position 4', 'Position 5'};
ndetectors = 6;
nbins = 500;

mu = plotEnergyFluence(sourcepath, filename, str, ndetectors, nbins)
