clear all;
close all;
fclose('all');
clc;

sourcefile = {'C:\Users\Delmon\Dropbox\PhD\FLUKA Monte Carlo\FLUKA files - 16 MeV proton dot grid\Dose.xlsx'};
str = {'16 MeV protons'};

for i = 1
    plotDose(sourcefile{i}, str{i})
end