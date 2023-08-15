clear all;
close all;
fclose('all');
clc;

%% Path
sourcepath_beam = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT in vitro\DCPT_beam_model.csv';
sourcepath_dcm  = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT in vitro\single1Gy_plan.dcm';
sourcepath_dose = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT in vitro\DCPT in vitro Dose LET\DCPT_invitro_21_plot.dat';
sourcepath_let  = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT in vitro\DCPT in vitro Dose LET\tab';
destpath        = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT in vitro';

%% Parameters
E0          = 0.106483;     % maximum energy in GeV % 0.027

lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y', ...
    '--r', '--g', '--c', '--m', '--k', '--y', ...
    ':r', ':g', ':c', ':m', ':k', ':y', ...
    '-.r', '-.g', '-.c', '-.m', '-.k', '-.y'};

%% Initialization
cnt                     = 0;
e                       = [];
N                       = [];
ScanSpotPositionMap     = {};
ScanSpotMetersetWeights = {};

%% Read info from DICOM RTPLAN

try
    % If dicominfo is successful, store the header information
    info = dicominfo(sourcepath_dcm, 'UseDictionaryVR', true);
catch
    % Otherwise, the file is either corrupt or not a real DICOM
    % file, so throw an error
    warning('File is not a valid DICOM object.');
end

for item = flip(fieldnames(info.IonBeamSequence.Item_1. ...
        IonControlPointSequence))'

    if any(strcmp(item{1}, {'Item_1', 'Item_3', 'Item_5', 'Item_7', ...
            'Item_9', 'Item_11', 'Item_13', 'Item_15'}))
        
        % Increment counter
        cnt = cnt + 1;

        % Store energy layers
        e(cnt) = info.IonBeamSequence.Item_1.IonControlPointSequence. ...
            (item{1}).NominalBeamEnergy/1000; % in GeV

        % Store number of scan spot positions
        N(cnt) = info.IonBeamSequence.Item_1.IonControlPointSequence. ...
            (item{1}).NumberOfScanSpotPositions;

        % Store positions (x-y pair coordinates) of the scan spots
        ScanSpotPositionMap{cnt} = reshape(double(info.IonBeamSequence. ...
            Item_1.IonControlPointSequence.(item{1}). ...
            ScanSpotPositionMap)./10, 2, [])';

        % Store Meterset weights corresponding to scan spot positions
        ScanSpotMetersetWeights_vec = double(info.IonBeamSequence.Item_1. ...
            IonControlPointSequence.(item{1}).ScanSpotMetersetWeights);
        CumulativeMetersetWeight(cnt) = sum(ScanSpotMetersetWeights_vec);
        ScanSpotMetersetWeights{cnt} = ScanSpotMetersetWeights_vec; % ./ ...
%             CumulativeMetersetWeight;

        % Find position and weight of central axis
        [~, idx] = pdist2([ScanSpotPositionMap{cnt}(:,1) ...
            ScanSpotPositionMap{cnt}(:,2)], [0 0], ...
            'euclidean', 'Smallest', 1);
        MetersetWeight_central(cnt) = ScanSpotMetersetWeights{cnt}(idx,:);

        % Plot scan spot configuration with according weights
        h = figure(cnt*10);
        set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        hold on
        plot(ScanSpotPositionMap{cnt}(:,2), ...
            ScanSpotPositionMap{cnt}(:,1), 'o--', 'LineWidth', 1, ...
            'MarkerSize', 5, 'MarkerEdgeColor', 'r', ...
            'MarkerFaceColor', [0.8500 0.3250 0.0980])
%         plot([-5 -5 5 5 -5], [-1.25 1.25 1.25 -1.25 -1.25], '-', ...
%             'LineWidth', 1, 'Color', 'k')
        hold off
%         xlim([-6.5 6.5])
%         ylim([-1.5 2.6])
        xlabel('x (cm)', 'Interpreter', 'latex', 'FontWeight','bold')
        ylabel('y (cm)', 'Interpreter', 'latex', 'FontWeight','bold')
        lgd = legend(['Energy layer $E_{nominal}=' num2str(e(cnt)*1000) '$ MeV'], ...
            'Location', 'NorthWest', 'Interpreter', 'latex');
        title(lgd, '\bf{Scan spot map w/ Meterset weights}')
        text(ScanSpotPositionMap{cnt}(:,2) + 0.1, ...
            ScanSpotPositionMap{cnt}(:,1) + 0.1, ...
            cellstr(num2str(round(ScanSpotMetersetWeights{cnt},4))), ...
            'FontSize', 9);
        set(gca, 'FontSize', 14)
        saveas(h, fullfile(destpath, ...
            sprintf('ScanSpotMap_E%iMeV.png', e(cnt)*10^3)));

    end

end

w = MetersetWeight_central ./ sum(MetersetWeight_central);

%% Get Gaussian parameters for the beam energies
opts = detectImportOptions(sourcepath_beam);
opts.SelectedVariableNames = {'x_EnergyNominal_MeV_', 'x_EReal_MeV_', ...
    'x_ERealSigma_MeV_', 'x_SigmaX_mm_', 'x_SigmaY_mm_', ...
    'x_SigmaX__rad_', 'x_SigmaY__rad_'};
T                       = readtable(sourcepath_beam, opts);
T.x_EnergyNominal_MeV_  = T.x_EnergyNominal_MeV_./1000; % in GeV
T.x_EReal_MeV_          = T.x_EReal_MeV_./1000;         % in GeV
T.x_ERealSigma_MeV_     = T.x_ERealSigma_MeV_./1000;    % in GeV
T.x_SigmaX_mm_          = T.x_SigmaX_mm_./10;           % in cm
T.x_SigmaY_mm_          = T.x_SigmaY_mm_./10;           % in cm
T.x_SigmaX__rad_        = T.x_SigmaX__rad_.*1000;       % in mrad
T.x_SigmaY__rad_        = T.x_SigmaY__rad_.*1000;       % in mrad

% Look-up table for 'real' energy
e_real = interp1(T.x_EnergyNominal_MeV_, T.x_EReal_MeV_, e', 'pchip')';

% Look-up table for sigma_e
sigma_e = interp1(T.x_EnergyNominal_MeV_, T.x_ERealSigma_MeV_, e', 'pchip')';
FWHM_e  = (2*sqrt(2*log(2))).*sigma_e;

% Look-up table for sigma_x_cm and sigma_y_cm
sigma_x_cm  = interp1(T.x_EnergyNominal_MeV_, T.x_SigmaY_mm_, e', 'pchip')';
sigma_y_cm  = interp1(T.x_EnergyNominal_MeV_, T.x_SigmaX_mm_, e', 'pchip')';
FWHM_x_cm   = (2*sqrt(2*log(2))).*sigma_x_cm;
FWHM_y_cm   = (2*sqrt(2*log(2))).*sigma_y_cm;

% Look-up table for sigma_x_mrad and sigma_y_mrad
sigma_x_mrad = interp1(T.x_EnergyNominal_MeV_, T.x_SigmaY__rad_, e', 'pchip')';
sigma_y_mrad = interp1(T.x_EnergyNominal_MeV_, T.x_SigmaX__rad_, e', 'pchip')';
FWHM_x_mrad  = (2*sqrt(2*log(2))).*sigma_x_mrad;
FWHM_y_mrad  = (2*sqrt(2*log(2))).*sigma_y_mrad;

% Use 'real' instead of nominal beam energy
e = e_real;

% [e'.*1000 sigma_e.*1000 sigma_x.*10 sigma_y.*10]
% mean(sigma_e)*1000

%% Write to file
% fid = fopen(fullfile(destpath, ...
%     ['SOBP_khi' num2str(khi*100) '_p' num2str(p*100) '.txt']), 'wt');
% fprintf(fid, '%.6f  %.6f   %.6f\n', [e_min; e_max; w]);
% fclose(fid);

fid = fopen(fullfile(destpath, ['BP_E0' num2str(e(end)*10^6) '.dat']), 'wt');
fprintf(fid, ['E[GeV]\t\tsigma_E[GeV]\tw_norm\t' ...
    'fwhm_x[cm]\tfwhm_y[cm]\tfwhm_x[mrad]\tfwhm_y[mrad]\n']);
fprintf(fid, ['%.6f\t%.6f\t %.6f\t %.6f\t %.6f\t %.6f\t%.6f\n'], ...
    [e; sigma_e; w; FWHM_x_cm; FWHM_y_cm; FWHM_x_mrad; FWHM_y_mrad]);
fclose(fid);

titlestr = ['Pristine Bragg peak (BP); ' ...
    '$E_{real}=' num2str(round(e(end)*1000,1)) '$ MeV, ' ...
    '$\sigma_E=' num2str(round(sigma_e(end)*1000,1)) '$ MeV, ' ...
    '$\sigma_x=' num2str(round(sigma_x_cm(end),2)) '$ cm, ' ...
    '$\sigma_y=' num2str(round(sigma_y_cm(end),2)) '$ cm'];

%% Averaged LET measurement vs depth

n_detectors     = length(getAllFiles(sourcepath_let)); % number of detectors (LET measurements)
n_bins          = 1000;                 % number of bins
depth_let       = linspace(1, 9, 70);   % depth in cm
depth_let       = [0 depth_let 9.1159];   

% Estimate averaged let
[LET_f, LET_d] = averageLET(sourcepath_let);

% Interpolate
depth_let_interp    = depth_let(1):0.005:depth_let(end);
LET_f               = interp1(depth_let, LET_f, depth_let_interp, 'pchip');
LET_d               = interp1(depth_let, LET_d, depth_let_interp, 'pchip');
depth_let           = depth_let_interp;

% Index of relevant depth interval
region  = 0 <= depth_let & depth_let <= 8.4; 
idx     = ismember(depth_let, [2 2.5 3 3.5 4 4.5 5 6.75 7.25 7.75 8.25]);

%% Read depth dose data

[depth_dose, dose] = getDepthDose(sourcepath_dose);

% Normalize to max
dose = (dose ./ max(dose(:))) * 100;

%% Plot
figure()
yyaxis left
plot(depth_dose, dose, 'k-', 'LineWidth', 1.0)
xlabel('Depth (cm)')
ylabel('Dose (%)') % nGy per primary particles/cm^2 per primary
xlim([0 inf])
ylim([0 inf])

yyaxis right
hold on
plot(depth_let(region), LET_f(region), 'b-', 'LineWidth', 1.0)
plot(depth_let(region), LET_d(region), 'r-', 'LineWidth', 1.0)
plot(depth_let(idx), LET_d(idx), 'ro')
hold off
ylabel('Mean LET (keV/\mum)')
ylim([0 inf])

legend('Dose', 'Fluence-averaged LET, $LET_f$', ....
    'Dose-averaged LET, $LET_d$', '$LET_d$ measurements', ...
    'Interpreter', 'latex', 'FontWeight', 'bold', 'Location', 'northwest')
title(titlestr, 'Interpreter', 'latex', 'FontWeight', 'bold')

grid on

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', 14)


% %% Plot LET spectrum at a selected depth
% figure()
% hold on
% for i = [15 63 65 67 68] %  1:n_detectors % 27
%     
%     if ~any(isnan(f(:,i)))     
%         plot(LETpoint(:,i), LETpoint(:,i) .* f(:,i))
%     end
%     
% end
% xlabel('L (keV/\mum)')
% ylabel('L * f(L) (keV/\mum)')
% % legend(['depth = ' num2str(depth(1) + (i-1)*median(diff(depth))) ' cm'])
% legend('Position 1', 'Position 2', 'Position 3', 'Position 4', 'Position 5')
