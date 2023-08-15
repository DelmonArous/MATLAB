clear all;
close all;
fclose('all');
clc;

sourcepath = {'C:\Users\Delmon\Dropbox\PhD\FLUKA Monte Carlo\FLUKA files - 16 MeV proton dot grid\USRYIELD - X transversal', ...
    'C:\Users\Delmon\Dropbox\PhD\FLUKA Monte Carlo\FLUKA files - 16 MeV proton dot grid\USRYIELD', ...
    'C:\Users\delmo\Desktop\AM paper\LET\tab'};
destpath = {'C:\Users\Delmon\Dropbox\PhD\FLUKA Monte Carlo\FLUKA files - 16 MeV proton dot grid\USRYIELD - X transversal\LET_16MeVproton.xlsx', ...
    'C:\Users\Delmon\Dropbox\PhD\FLUKA Monte Carlo\FLUKA files - 16 MeV proton dot grid\USRYIELD\LET_16MeVproton.xlsx', ...
    'C:\Users\delmo\Desktop\AM paper'};

%% Data sheet structure
n_detectors = [73 70 70];  % number of detectors (LET measurements)
n_bins = 1000;          % number of bins

%% Depths of each averaged LET measurement

% 16 MeV proton
% depth = 0:0.005:0.25;
% add1 = 0.15125:0.005:0.19625;
% add2 = 0.1525:0.005:0.1975;
% add3 = 0.15375:0.005:0.19875;
% depth = sort([depth add1 add2 add3 0.0008]);
%depth = sort([depth 0.0008]);

% 16 MeV proton Falcon
% depth = 0:0.005:0.205;   
% add1 = 0.15125:0.005:0.19625;
% add2 = 0.1525:0.005:0.1975;
% add3 = 0.15375:0.005:0.19875;
% depth = sort([depth add1 add2 add3]);

% depth = 2.6:0.0025:2.7725;
depth = linspace(81.287-(81.3-0.013), 81.4-(81.3-0.013), 70);

%depth = -1.8:0.0125:-0.9;

%%
for i = 3
    LETaverageEstimate(sourcepath{i}, destpath{i}, n_detectors(i), n_bins, depth)
end