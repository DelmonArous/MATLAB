clear all;
close all;
clc;

% Read in data-file
data = load('dlts_A_4E10_down.txt');
T = data(:,1);
C_r = data(:,2); % A-4e10 down
w1 = data(:,4);
w2 = data(:,5);
w3 = data(:,6);
w4 = data(:,7);
w5 = data(:,8);
w6 = data(:,9);

% C_r = data(:,4);
% w1 = data(:,6);
% w2 = data(:,7);
% w3 = data(:,8);
% w4 = data(:,9);
% w5 = data(:,10);
% w6 = data(:,11);

% Plot of the DLTS-spectra for time window 1-6
figure()
plot(T, w1./C_r, '-r', T, w2./C_r, '-g', T, w3./C_r, '-b',...
    T, w4./C_r, '-c', T, w5./C_r, '-m', T, w6./C_r, '-k')
legend('Window 1', 'Window 2', 'Window 3',...
    'Window 4', 'Window 5', 'Window 6');
xlabel('Temperature T [K]','interpreter','latex','FontSize',13)
ylabel(['$DLTS/C_r$'],'interpreter','latex','FontSize',13)
% Plot of the DLTS-spectrum for time window 6
% figure()
% plot(T, DLTS_w6./0.102)
% legend('Window 6')
% xlabel(['$Temperature\,[K]$'],'interpreter','latex','FontSize',13)
% ylabel(['$DLTS/0.102C_r\,$'],'interpreter','latex','FontSize',13)

% Important constants
N_d = 3.0*10^13; % doping concentration [cm^-3]
% Emission rate constants for the different time windows [w1 w2 ... w6]
e_n = [70.66836 38.64953 20.28329 10.40057 5.267648 2.651011];
k_J = 1.3806488*10^(-23); % Boltzmann constant [J/K]
k_eV = 8.6173324*10^(-5); % Boltzmann constant [eV/K]
h = 6.62606957*10^(-34); % Planck constant [Js]
m_e = 9.10938291*10^(-31); % electron mass [kg]
beta = sqrt((3*k_J)/(1.06*m_e))*2*((2*pi*1.06*m_e*k_J)/h^2)^(3/2)*10^(-4);

for i = 1:3 % loops over each defect (peak)
   if (i == 1) % 1.peak
       T_peak = [90. 88. 87. 84.  82. 80.]; % A - 4e10
       %T_peak = [93. 91. 90. 89. 87. 84.]; % A - 2e10
       %T_peak = [92. 91. 89. 87. 86. 84.]; % B - 4e10
       %T_peak = [92. 89. 88. 86. 84. 83.]; % B - 3e9
       %T_peak = [91. 90. 88. 87. 86. 84.]; % B - 2e10
       %T_peak = [94. 92. 90. 88.5 87. 85.];
       
       % S_max read from sample [B_3e9_up B_6e9_up B_2e10_up A_2e10]
       S_max = [0.006003 0.01326 0.04654 0.04681];
       dose = [3e9 6e9 2e10 2e10]; % irradiated proton doses
   elseif (i == 2) % 2.peak 
       T_peak = [129. 125. 122. 120. 117. 114.]; % A - 4e10
       %T_peak = [130. 128. 123. 120. 118. 115.]; % A - 2e10
       %T_peak = [132. 128. 125. 121. 121. 116.]; % B - 4e10
       %T_peak = [130. 127. 129. 122. 118. 115.]; % B - 3e9
       %T_peak = [128. 126. 124. 121. 117. 114.]; % B - 2e10
       %T_peak = [123. 119. 118. 116. 112. 111.]; % S1 - 6e9
       %T_peak = [130. 127. 124. 121. 118. 115.];
       
       % S_max read from sample [B_3e9_up B_6e9_up B_2e10_up A_2e10 B_4e10_up A_4e10_down]
       S_max = [0.002244 0.004044 0.01295 0.01293 0.02001 0.01728];  
       dose = [3e9 6e9 2e10 2e10 4e10 4e10];
   elseif (i == 3) % 3.peak
       T_peak = [222. 215. 211. 207. 200. 196.]; % A - 4e10
       %T_peak = [217. 216. 212. 205. 203. 196.]; % A - 2e10
       %T_peak = [220. 215. 214. 208. 204. 199.]; % B - 4e10
       %T_peak = [222. 220. 213. 207. 203. 197.]; % B - 3e9
       %T_peak = [223. 219. 213. 210. 204. 199.5]; % B - 2e10
       %T_peak = [220. 222. 217. 209. 205. 200.]; % S1 - 6e9
       
       % S_max read from sample [B_3e9_up B_6e9_up B_2e10_up A_2e10 B_4e10_up A_4e10_down]
       S_max = [0.002575 0.005152 0.01698 0.01818 0.03292 0.03539];
       dose = [3e9 6e9 2e10 2e10 4e10 4e10];
   end
   
   % Arrhenius plot
   figure()
   plot(1./T_peak, log(e_n./T_peak.^2), '*')
   xlabel(['$1/T\,[K^{-1}]$'],'interpreter','latex','FontSize',13)
   ylabel('$\ln\left(e_n/T^2\right)$','interpreter','latex','FontSize',13)
   p1 = polyfit(1./T_peak, log(e_n./T_peak.^2), 1);
   yfit1 = polyval(p1, 1./T_peak);
   hold on
   plot(1./T_peak, yfit1, '-ro')
   str1 = sprintf('Linefit: %0.5g + %0.5g T^{-2}', p1(2), p1(1));
   legend('Data', str1)
   
   % Calculation of the electron state DeltaE =  E_c - E_T
   E_T = -p1(1)*8.6173325*10^(-5)
   % Calculation of the trap cross-section sigma_n(p)
   sigma = exp(p1(2))/beta;
   
   % Calculation of the trap concentration N_T 
   % and plot of N_T vs irradiated proton dose 
   figure()
   N_T = N_d*S_max;
   N_T = N_T./N_d;
   plot(dose, N_T, '*')
   xlabel(['Dose $[H^+/cm^2]$'],'interpreter','latex','FontSize',13)
   ylabel(['$N_T\,[cm^{-3}]$'],'interpreter','latex','FontSize',13)
   p2 = polyfit(dose, N_T, 1);
   yfit2 = polyval(p2, dose);
   hold on
   plot(dose, yfit2, '-ro')
   str2 = sprintf('Linefit: %0.5g + %0.5g x_{dose}', p2(2), p2(1));
   legend('Data', str2);
end