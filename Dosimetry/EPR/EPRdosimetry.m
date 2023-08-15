clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')


%% Path

filename = 'C:\Users\delmo\Desktop\EPR dosimetry\EPR_dosimetry_160921.xlsx';

%% Read calibration data

for i = 3 % 3
   
    data_calib    = xlsread(filename, i);
    dose_calib    = data_calib(:,1);
    PP_calib      = data_calib(:,end);
    
end

%% Read measurement data

data_Open   = xlsread(filename, 4);
data_GRID = xlsread(filename, 5);

PP_Open     = data_Open(:,end);
PP_GRID     = data_GRID(:,end);

position_Open = {'Control'; 'Control'; 'Control'; 'Control'; 'Control'; ...
    'A'; 'A'; 'A'; 'A'; 'A'; 'B'; 'B'; 'B'; 'B'; 'B'; 'C'; 'C'; 'C'; ...
    'C'; 'C'; 'D'; 'D'; 'D'; 'D'; 'D'};
position_GRID = {'A_P'; 'A_P'; 'A_P'; 'A_P'; 'A_V'; 'A_V'; 'A_V'; 'A_V'; ...
    'A_V'; 'A_V'; 'B_P'; 'B_P'; 'B_P'; 'B_P'; 'B_V'; 'B_V'; 'B_V'; ...
    'B_V'; 'C_P'; 'C_P'; 'C_P'; 'C_P'; 'C_V'; 'C_V'; 'C_V'; 'C_V'; ...
    'D_P'; 'D_P'; 'D_P'; 'D_P'; 'D_V'; 'D_V'; 'D_V'; 'D_V'};

%% Simple linear regression

% X = [ones(length(dose_calib),1) dose_calib];
% b = X\PP_calib;
% 
% fit = X*b;

tbl         = table(dose_calib, PP_calib, 'VariableNames', {'Dose','PP'});
lm          = fitlm(tbl, 'PP~Dose');

% Linear regression fit estimate
fit_calib   = lm.Coefficients.Estimate(1) + ...
    lm.Coefficients.Estimate(2).* dose_calib;

% Estimate doses
dose_Open = (PP_Open - lm.Coefficients.Estimate(1)) ./ ...
    lm.Coefficients.Estimate(2);
dose_GRID = (PP_GRID - lm.Coefficients.Estimate(1)) ./ ...
    lm.Coefficients.Estimate(2);

tbl_Open = table(dose_Open, position_Open, ...
    'VariableNames', {'Dose','Position'})
tbl_GRID = table(dose_GRID, position_GRID, ...
    'VariableNames', {'Dose','Position'})

peakind     = [1:4 11:14 19:22 27:30];
valleyind   = [5:10 15:18 23:26 31:34];
[mean(dose_Open(6:end)) std(dose_Open(6:end))]
[mean(dose_GRID(peakind)) std(dose_GRID(peakind))]
[mean(dose_GRID(valleyind)) std(dose_GRID(valleyind))]

%% Plot

figure();
plot(dose_calib, PP_calib, 'o', dose_calib, fit_calib, '--')
legend('Data', 'Linear fit', 'Location', 'NorthWest');
xlabel('Dose (Gy)')
ylabel('Peak-to-Peak')
xlim([min(dose_calib) max(dose_calib)+1])
grid on
set(gca, 'FontSize', 16)
