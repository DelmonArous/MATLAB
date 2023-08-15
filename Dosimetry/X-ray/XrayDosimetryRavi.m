clc
clear all
close all

%% Date

date = datestr(datetime(floor(now), 'ConvertFrom', 'datenum'));
date = strrep(date, '-', '_');

%% Beam parameters 
% E       = 225; % in kV
% current = 10; % in mA 
% filter  = '0.7 mm Cu + 1.52 mm Al';
% SSD     = 50; % in cm

%% Parameters
N_k     = 43.77/1000;    % in Gy/nC
SD_N_k  = 0.39/1000;     % in Gy/nC
K_u     = 1.0;
mu_en   = 1.075;
P_u     = 1.02;
T0      = 20.0;     % in degrees Celcius
P0      = 1013;     % in hPa
T       = 22.2;     % in degrees Celcius
P       = 1025;      % in hPa
k_TP    = ((273.15 + T)/ (273.15 + T0)) * (P0/P);

dose_nom = [0.3; 0.5; 0.6; 0.9; 1.0; 1.2; 1.5; 1.8];   % nominal dose in Gy
t_rad = [5; 10; 15; 20]; % low dosimetry irradiation time in sec
D_high_vec      = [];
t_estimate_vec  = [];

%% M_u in nC
M_u_low     = {[0.40; 0.48; 0.65; 0.56], [1.81; 2.03; 1.96; 1.79], ...
    [3.45; 3.20; 3.34; 3.28], [4.67; 4.77; 4.79; 4.65]}; % 5, 10, 15, 20 sec; 4 measurements each
M_u_high    = [15.71; 15.88; 15.90; 15.55];   % 60 sec; 4 measurements each
% save(['M_u_high' date '_new.mat'], 'M_u_high')
% save(['M_u_low' date '_new.mat'], 'M_u_low')

%%

% Low dose range dosimetry
b = LowDoseXrayDosimetry(M_u_low, N_k, K_u, mu_en, P_u, k_TP, SD_N_k, t_rad);
intercept   = b(1);  
slope       = b(2);

% High dose range dosimetry
[D_w_high, SD_D_w_high] = ...
    DoseToWater(M_u_high, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
D_high_vec = [D_high_vec; D_w_high];

%% Estimate irradiation time

% D_high_vec
D_high          = mean(D_high_vec); 

t_estimate = ((dose_nom(1:3) - intercept)./slope)/60;
t_estimate_vec = [t_estimate_vec; t_estimate];
t_estimate      = dose_nom(4:end)./D_high;
t_estimate_vec  = [t_estimate_vec; t_estimate];
t_estimate_vec  = minutes(t_estimate_vec);
t_estimate_vec.Format = 'mm:ss'

