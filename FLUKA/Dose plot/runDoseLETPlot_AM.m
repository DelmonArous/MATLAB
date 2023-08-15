clear all;
close all;
fclose('all');
clc;

%% Directories 
sourcepath_dose = ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Manuscripts\AM paper\Dose\Proton 15.22 MeV 103.1 cm Dist\Dose.dat';
sourcepath_LET = ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Manuscripts\AM paper\LET\Proton 15.22 MeV 103.1 cm Dist\tab';

%% Variables

n_detectors = 70;   % number of detectors (LET measurements)
n_bins = 1000;      % number of bins

% lgdtitle = {'Lid, 81.3 cm', 'Parafilm, 81.3 cm'};
% str = {'Dose', 'Measured dose', 'LET_d', 'LET_f'};
str = {'Dose', 'LET_d', 'Measurements'};
%str = {'Dose', 'LET_d', 'LET_f'};
lgdtitle = {'15.22 MeV proton'};

lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y', ...
    '--r', '--g', '--c', '--m', '--k', '--y', ...
    ':r', ':g', ':c', ':m', ':k', ':y', ...
    '-.r', '-.g', '-.c', '-.m', '-.k', '-.y'};

corr_dist = 81.3-0.013;
depth_LET = linspace(81.287-corr_dist, 81.4-corr_dist, n_detectors); % 81.287

%% Read depth dose data

[filepath, filename, ext] = fileparts(sourcepath_dose);
[depth_dose, Dose] = plotDose(fullfile(filepath, [filename ext]));

% Convert to Gy
%         Dose = Dose .* 1.602176462*10^(-7) * 10^9; % in nGy
% Normalize to max
Dose = (Dose ./ max(Dose(:))) * 100;
% Correct file depth
depth_dose = depth_dose - corr_dist; % corr_dist(j);

%% Read LET data

[LET_d, LET_f] = LETaverageEstimate(sourcepath_LET, n_detectors, depth_LET);

% [LET_d_test, LET_f_test] = LETaverageEstimate('C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\AM paper\LET\Proton 15.22 MeV\v2', ...
%     5, depth_LET);
% 
% LET_d_test 
% % LET_f_test
% 
% [LET_d_test, LET_f_test] = LETaverageEstimate('C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\AM paper\LET\Proton 15.22 MeV\v3', ...
%     5, depth_LET);
% 
% LET_d_test 
% % LET_f_test
% 
% [LET_d_test, LET_f_test] = LETaverageEstimate('C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\AM paper\LET\Proton 15.22 MeV\v4', ...
%     5, depth_LET);

% LET_d_test 
% LET_f_test

%% Plot

% [0.0162 0.1132] % [Pos1 Pos5] 
% [0.0133 0.1006 0.1041 0.1069 0.1091]

figure()
yyaxis left
% plot(depth_dose, Dose, lineSpec{1}, ...
%     [0.0133 0.1006 0.1041 0.1069 0.1091], ...
%     [23.5003056 100.0 82.0255921 50.319236 24.2433081], ...
%     'ro', 'LineWidth', 1.0)
plot(depth_dose, Dose, lineSpec{1}, 'LineWidth', 1.0)
% hold on
% errorbar([0.0236 0.1096], [23.5003056 24.2433081], ...
%     [0.144935435 0.431778107], 'vertical', 'ro')

xlabel('Depth (cm)')
ylabel('Dose (%)') % nGy per primary, particles/cm^2 per primary
xlim([0 inf])
ylim([0 inf])

yyaxis right
% plot(depth_LET, LET_d, lineSpec{2},  ...
%     depth_LET, LET_f, lineSpec{3}, 'LineWidth', 1.0)
plot(depth_LET, LET_d, lineSpec{3}, ...
    [0.02 0.102 0.1055 0.1075 0.11], ...
    [4.7 25.1 33.0 37.0 40.4], 'bo', 'LineWidth', 1.0)
ylabel('Mean LET (keV/\mum)')
ylim([0 inf])

lgd = legend(str, 'Location', 'NorthWest');
title(lgd, lgdtitle)

% plot(depth, LET_d, 'b-o', depth, LET_f, 'r-x') %, ...
%     %depth_interp, LET_d_interp, 'b-', depth_interp, LET_f_interp, 'r-')
% legend('LET_d', 'LET_f', 'Location', 'best')

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', 14)
