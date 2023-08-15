clc
clear all
close all

N_k     = 43.77/1000;    % in Gy/nC
SD_N_k  = 0.39/1000;     % in Gy/nC
K_u     = 1.0;
mu_en   = 1.075;
P_u     = 1.02;
T0      = 20.0;             % in degrees Celcius
P0      = 1013;             % in hPa
T       = 23.3; % 23.8;     % in degrees Celcius
P       = 1016; % 1004;      % in hPa
k_TP    = ((273.15 + T)/ (273.15 + T0)) * (P0/P);

dose_nom = [0.3; 0.6; 0.9; 1.2; 1.5; 1.8];   % nominal dose in Gy
t_rad = [5; 10; 15; 20];      % low dosimetry irradiation time in sec
intercept_vec   = [];
slope_vec       = [];
D_high_vec      = [];
SD_D_high_vec   = [];
t_estimate_vec  = [];

%% Pos A: 
M_u_low = {[0.49 0.33 0.38 0.43], ...   % 5 sec; 4 measurements each
            [1.53 1.53 1.46 1.48], ...  % 10 sec; 4 measurements each
            [2.58 2.42 2.43 2.56], ...  % 15 sec; 4 measurements each
            [3.54 3.57 3.45 3.51]};     % 20 sec; 4 measurements each
M_u_high = [12.04 12.00 12.05 12.14];   % 60 sec; 4 measurements each

% Low dose range dosimetry
b = LowDoseXrayDosimetry(M_u_low, N_k, K_u, mu_en, P_u, k_TP, SD_N_k, t_rad);
intercept_vec   = [intercept_vec; b(1)];  
slope_vec       = [slope_vec; b(2)];

% High dose range dosimetry
[D_w_high, SD_D_w_high] = ...
    DoseToWater(M_u_high, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
D_high_vec = [D_high_vec; D_w_high];
SD_D_high_vec = [SD_D_high_vec; SD_D_w_high];

%% Pos B:
M_u_low = {[0.42 0.36 0.47 0.42], ...   % 5 sec; 4 measurements each
            [1.44 1.46 1.39 1.49], ...  % 10 sec; 4 measurements each
            [2.51 2.33 2.49 2.40], ...  % 15 sec; 4 measurements each
            [3.52 3.42 3.47 3.38]};     % 20 sec; 4 measurements each
M_u_high = [12.25 12.33 12.13 12.29];   % 60 sec; 4 measurements each

% Low dose range dosimetry
b = LowDoseXrayDosimetry(M_u_low, N_k, K_u, mu_en, P_u, k_TP, SD_N_k, t_rad);
intercept_vec   = [intercept_vec; b(1)];  
slope_vec       = [slope_vec; b(2)];

% High dose range dosimetry
[D_w_high, SD_D_w_high] = ...
    DoseToWater(M_u_high, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
D_high_vec = [D_high_vec; D_w_high];
SD_D_high_vec = [SD_D_high_vec; SD_D_w_high];

%% Pos C:
M_u_low = {[0.47 0.33 0.31 0.44], ...   % 5 sec; 4 measurements each
            [1.61 1.52 1.53 1.52], ...  % 10 sec; 4 measurements each
            [2.51 2.60 2.57 2.62], ...  % 15 sec; 4 measurements each
            [3.61 3.64 3.57 3.68]};     % 20 sec; 4 measurements each
M_u_high = [11.95 11.86 11.91 11.94];   % 60 sec; 4 measurements each

% Low dose range dosimetry
b = LowDoseXrayDosimetry(M_u_low, N_k, K_u, mu_en, P_u, k_TP, SD_N_k, t_rad);
intercept_vec   = [intercept_vec; b(1)];  
slope_vec       = [slope_vec; b(2)];

% High dose range dosimetry
[D_w_high, SD_D_w_high] = ...
    DoseToWater(M_u_high, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
D_high_vec = [D_high_vec; D_w_high];
SD_D_high_vec = [SD_D_high_vec; SD_D_w_high];

%% Pos D:
M_u_low = {[0.41 0.30 0.37 0.43], ...   % 5 sec; 4 measurements each
            [1.47 1.56 1.45 1.37], ...  % 10 sec; 4 measurements each
            [2.56 2.61 2.47 2.40], ...  % 15 sec; 4 measurements each
            [3.56 3.52 3.50 3.68]};     % 20 sec; 4 measurements each
M_u_high = [11.96 11.92 12.10 11.94];   % 60 sec; 3 measurements each

% Low dose range dosimetry
b = LowDoseXrayDosimetry(M_u_low, N_k, K_u, mu_en, P_u, k_TP, SD_N_k, t_rad);
intercept_vec   = [intercept_vec; b(1)];  
slope_vec       = [slope_vec; b(2)];

% High dose range dosimetry
[D_w_high, SD_D_w_high] = ...
    DoseToWater(M_u_high, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
D_high_vec = [D_high_vec; D_w_high];
SD_D_high_vec = [SD_D_high_vec; SD_D_w_high];

%% Estimate irradiation time

intercept_vec
slope_vec
D_high_vec
SD_D_high_vec

intercept   = mean(intercept_vec);
slope       = mean(slope_vec);
D_high      = mean(D_high_vec); 


% Scaling factor: 1+(1-0.9744)
t_estimate = ((dose_nom(1).*1.0256 - intercept)./slope)/60;
t_estimate_vec = [t_estimate_vec; t_estimate];
t_estimate = (dose_nom(2:end).*1.0256)./D_high;
t_estimate_vec = [t_estimate_vec; t_estimate];

t_estimate_vec = minutes(t_estimate_vec);
t_estimate_vec.Format = 'mm:ss'
