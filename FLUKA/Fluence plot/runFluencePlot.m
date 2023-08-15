clear all;
close all;
fclose('all');
clc;

sourcepath = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\X-ray 220 keV\Fluence';
str = {{'Secondary electrons, 1.52 mm Al + 0.70 mm Cu filter', ...
    'Secondary electrons, no filter', ...
    'Primary photons, 1.52 mm Al + 0.70 mm Cu filter', ...
    'Primary photons, no filter'}};

plotFluence(sourcepath, str{1})