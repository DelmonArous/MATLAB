clear all;
close all;
fclose('all');
clc;

% sourcepath = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\X-ray 220 keV\Dose';
% str = {{'220 kV X-ray, 1.52 mm Al + 0.70 mm Cu filter', ...
%     '220 kV X-ray, no filter'}};

sourcepath = 'C:\Users\delmo\Desktop\AM paper\Dose\Proton 15.2 MeV Water target\XLSX';
str = {{'15.2 MeV proton'}};

plotDose_v2(sourcepath, str{1})