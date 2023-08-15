clc
clear all
close all

N_k     = 43.77/1000;    % in Gy/nC
SD_N_k  = 0.39/1000;     % in Gy/nC
K_u     = 1.0;
mu_en   = [1.027 1.027 1.027 1.068 1.068 1.068 1.075 1.075 1.075 1.075 ...
    1.075 1.027 1.027 1.027 1.068 1.068 1.068 1.027 1.027 1.068 1.068 ...
    1.027 1.027 1.027 1.027 1.027 1.027 1.027 1.027 1.027 1.027];
P_u     = [1.03 1.03 1.03 1.02 1.02 1.02 1.02 1.02 1.02 1.02 1.02 1.03 ...
    1.03 1.03 1.02 1.02 1.02 1.03 1.03 1.02 1.02 1.03 1.03 1.03 1.03 ...
    1.03 1.03 1.03 1.03 1.03 1.03];
T0      = 20.0;     % in degrees Celcius
P0      = 1013;     % in hPa
T       = 27.0;     % in degrees Celcius
P       = 1020-0.13;     % in hPa
k_TP    = ((273.15 + T)/ (273.15 + T0)) * (P0/P);

M_u = {[82.08 81.59], ...       100 kV, 2 mm Al, 10 mA, 5.00 Gy nominal dose
    [106.01 106.06], ...        100 kV, 2 mm Al, 10 mA, 6.50 Gy nominal dose
    [138.37 138.25], ...        100 kV, 2 mm Al, 10 mA, 8.50 Gy nominal dose
    [88.23 87.99 87.98], ...    180 kV, 0.3 mm Cu, 10 mA, 5.00 Gy nominal dose
    [101.17 101.00 100.87], ... 180 kV, 0.3 mm Cu, 10 mA, 5.75 Gy nominal dose
    [114.26 114.54 113.82], ... 180 kV, 0.3 mm Cu, 10 mA, 6.50 Gy nominal dose
    [172.38 172.42 172.99], ...  225 kV, 0.3 mm Cu, 10 mA, 8.50 Gy nominal dose
    [200.24 201.43 202.35 201.43 202.00], ... 225 kV, 0.3 mm Cu, 10 mA, 8.50 Gy nominal dose
    [185.58 177.98 185.70], ... 225 kV, 0.3 mm Cu, 10 mA, 8.50 Gy nominal dose
    [183.62 182.99 183.54], ... 225 kV, 0.3 mm Cu, 10 mA, 8.50 Gy nominal dose, Pos: in box, right skewed
    [181.31 182.01 183.07], ... 225 kV, 0.3 mm Cu, 10 mA, 8.50 Gy nominal dose, Pos: in box, right straight
    [83.17 83.23 83.07], ...    100 kV, 2.0 mm Al, 10 mA, 5.00 Gy nominal dose, Pos: in box
    [95.54 96.58 96.67], ...    100 kV, 2.0 mm Al, 10 mA, 5.75 Gy nominal dose, Pos: in box
    [108.56 109.01 109.15], ... 100 kV, 2.0 mm Al, 10 mA, 6.50 Gy nominal dose, Pos: in box
    [89.28 89.80 89.28], ...    180 kV, 0.3 mm Cu, 10 mA, 5.00 Gy nominal dose, Pos: in box
    [102.97 103.12 102.48], ... 180 kV, 0.3 mm Cu, 10 mA, 5.75 Gy nominal dose, Pos: in box
    [115.42 115.73 116.12], ... 180 kV, 0.3 mm Cu, 10 mA, 6.50 Gy nominal dose, Pos: in box
    [84.11 83.53], ...          100 kV, 2.0 mm Al, 15 mA, 5.00 Gy nominal dose, Pos: in box
    [109.43 109.28], ...        100 kV, 2.0 mm Al, 15 mA, 6.50 Gy nominal dose, Pos: in box
    [89.35 89.42], ...          180 kV, 0.3 mm Cu, 15 mA, 5.00 Gy nominal dose, Pos: in box
    [116.22 115.99], ...        180 kV, 0.3 mm Cu, 15 mA, 6.50 Gy nominal dose, Pos: in box
    [148.18 148.37 148.33], ... 100 kV, 2.0 mm Al, 15 mA, 8.50 Gy nominal dose, Pos: in box (right) 
    [147.40 142.86 142.61], ... 100 kV, 2.0 mm Al, 15 mA, 8.50 Gy nominal dose, Pos: in box (left)
    [13.57 13.57 13.52], ...    100 kV, 2.0 mm Al, 15 mA, 1 min irradiation (0.740, 0.739, 0.752 Gy nominal dose), Pos: in box (right)
    [13.55 13.58 13.61] ...     100 kV, 2.0 mm Al, 15 mA, 1 min irradiation (0.742, 0.753, 0.741 Gy nominal dose), Pos: in box (left)
    [129.56 127.42 125.59 127.96] ...           100 kV, 2.0 mm Al, 15 mA, 7.0 Gy nominal dose, Pos: in box (right)
    [123.40 123.56 123.70 122.79 129.64] ...    100 kV, 2.0 mm Al, 15 mA, 7.0 Gy nominal dose, Pos: in box (left)
    [36.37 36.36 36.33] ...                     100 kV, 2.0 mm Al, 15 mA, 2.0 Gy nominal dose, Pos 6 (right)
    [35.62 35.53 35.76] ...                     100 kV, 2.0 mm Al, 15 mA, 2.0 Gy nominal dose, Pos 6 (left)
    [34.07 34.06 33.50] ...                     100 kV, 2.0 mm Al, 15 mA, 2.0 Gy nominal dose, Pos 4 (right)
    [31.98 32.10 31.99] ...                     100 kV, 2.0 mm Al, 15 mA, 2.0 Gy nominal dose, Pos 4 (left)
    };

dose_nominal = [5.0 6.5 8.5 5.0 5.75 6.5 8.5 8.5 8.5 8.5 8.5 5.0 5.75 ...
    6.5 5.0 5.75 6.5 5.0 6.5 5.0 6.5 8.5 8.5 ...
    mean([0.740 0.739 0.752]) mean([0.742, 0.753, 0.741]) ...
    7.0 7.0 2.0 2.0 2.0 2.0];

for i = 26:length(M_u) 
    
    [D_w, SD_D_w] = ...
        DoseToWater(M_u{i}, N_k, K_u, mu_en(i), P_u(i), k_TP, SD_N_k);
    [mean(M_u{i}) std(M_u{i})]
    [D_w, SD_D_w]
    dose_nominal(i)./D_w
    
end


% % Position 6, 3.5 Gy nominal dose
% M_u_pos6    = [61.50 61.21 60.95 60.99];
% 
% % Position 5, 3.5 Gy nominal dose
% M_u_pos5    = [61.33 60.97 60.67];
% 
% % Position 4, 3.5 Gy nominal dose
% M_u_pos4    = [58.48 58.16 58.51];
% 
% % Position 6, 5.0 Gy nominal dose
% M_u_5Gy     = [88.08];  
% 
% % Calculate dose to water
% D_w6    = mean(M_u_pos6) * N_k * K_u * mu_en * P_u * k_TP;
% D_w5    = mean(M_u_pos5) * N_k * K_u * mu_en * P_u * k_TP;
% D_w4    = mean(M_u_pos4) * N_k * K_u * mu_en * P_u * k_TP;
% D_w5Gy  = mean(M_u_5Gy) * N_k * K_u * mu_en * P_u * k_TP;
% 
% % Calculate corresponding standard deviation 
% SD_D_w6 = D_w6 * sqrt((SD_N_k/N_k)^2 + (std(M_u_pos6)/mean(M_u_pos6))^2);
% SD_D_w5 = D_w5 * sqrt((SD_N_k/N_k)^2 + (std(M_u_pos5)/mean(M_u_pos5))^2);
% SD_D_w4 = D_w4 * sqrt((SD_N_k/N_k)^2 + (std(M_u_pos4)/mean(M_u_pos4))^2);
% SD_D_w5Gy = D_w5Gy * sqrt((SD_N_k/N_k)^2 +(std(M_u_5Gy)/mean(M_u_5Gy))^2);






% %% Jacob, striped GRID
% 
% N_k     = 43.77/1000;    % in Gy/nC
% SD_N_k  = 0.39/1000;     % in Gy/nC
% K_u     = 1.0;
% mu_en   = 1.075;
% P_u     = 1.02;
% T0      = 20.0;     % in degrees Celcius
% P0      = 1013.25;  % in hPa
% T       = 23.5;     % in degrees Celcius
% P       = 1018-0.13;     % in hPa
% k_TP    = ((273.15 + T)/ (273.15 + T0)) * (P0/P); 
% 
% % 220 kV, 10 mA, 0.7 mm Cu + 1.52 mm Al, 1 min irradiation
% M_u = { [11.82 11.68 11.71], ...    % Pos A
%         [11.70 11.65 11.61], ...    % Pos B
%         [11.88 11.86 11.92], ...    % Pos C
%         [12.11 12.02 12.08]         % Pos D
%         };
% 
% % dose_nom = [5.0 6.5 8.5 5.0 5.75 6.5 8.5 8.5 8.5 8.5 8.5 5.0 5.75 6.5 5.0 5.75 6.5 5.0 6.5 5.0 6.5];
% 
% for i = 1:length(M_u)
%     
%     [D_w, SD_D_w] = ...
%         DoseToWater(M_u{i}, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
%     [D_w/60, SD_D_w]
%     
%     %     [mean(M_u{i}) std(M_u{i})]
% %     dose_nom(i)/D_w
%     
% end

