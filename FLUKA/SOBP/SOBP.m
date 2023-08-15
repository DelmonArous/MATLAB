clear all;
close all;
fclose('all');
clc;

%% Path
sourcepath_beam = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT\DCPT_beam_model.csv';
sourcepath_SOBP = {'C:\Users\delmo\Desktop\SOBP\Jette\Dose', ...
    'C:\Users\delmo\Desktop\SOBP\Exp\Dose_v2', ...
    'C:\Users\delmo\Desktop\SOBP\Beamlets\Dose'};
sourcepath_dcm  = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA Monte Carlo\DCPT\plan_22Gy_phys.dcm';
destpath        = 'C:\Users\delmo\Desktop\SOBP';

%% Parameters
E0          = 0.106483;     % maximum energy in GeV % 0.027
sigma_E0    = 0.0011144;    % Gaussian energy spread in GeV
alpha       = 0.0022;
p0          = 1.77;
n           = 7;            % number of energy intervals (layers) in SOBP
chi         = 0.35;         % fraction of the full range taken by the SOBP
p           = 1.40;
% 1.35:0.01:1.40 loop over different values of p for khi=0.40
% Possibly increase khi to >0.40 to achieve homogeniety in scoring volume
% -> then loop over different values of p

lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y', ...
    '--r', '--g', '--c', '--m', '--k', '--y', ...
    ':r', ':g', ':c', ':m', ':k', ':y', ...
    '-.r', '-.g', '-.c', '-.m', '-.k', '-.y'};
titlestr = ['Spread-out Bragg peak (SOBP); $E_0=' num2str(E0*1000) '$ MeV, ' ...
    '$\chi=' num2str(chi*100) '\%$, $\alpha=' num2str(alpha*100) '$, ' ...
    '$p=' num2str(p) '$, $n=' num2str(n+1) '$ energy layers, ' ...
    'Jette D. \& Chen W., (2011)'];

%% Full range
R0 = R(E0, alpha, p0);

%% Initialization
r = [];
e = [];
w = [];

%% Loop through proton energies

for k = 0:n

    r(k+1) = (1 - (1 - (k/n)) * chi) * R0;
    e(k+1) = (r(k+1)/alpha)^(1/p0);

    if k == 0
        w(k+1) = 1 - (1 - (1/(2*n)))^(1 - (1/p));
    elseif k == n
        w(k+1) = (1/(2*n))^(1 - (1/p));
    else
        w(k+1) = (1 - (1/n)*(k - 0.5))^(1 - (1/p)) - ...
            (1 - (1/n)*(k + 0.5))^(1 - (1/p));
    end

end

%% Read info from DICOM RTPLAN

try
    % If dicominfo is successful, store the header information
    info = dicominfo(sourcepath_dcm, 'UseDictionaryVR', true);
catch
    % Otherwise, the file is either corrupt or not a real DICOM
    % file, so throw an error
    warning('File is not a valid DICOM object.');
end

% Initialize
cnt                     = 0;
ScanSpotPositionMap     = {};
ScanSpotMetersetWeights = {};

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
        plot([-5 -5 5 5 -5], [-1.25 1.25 1.25 -1.25 -1.25], '-', ...
            'LineWidth', 1, 'Color', 'k')
        hold off
        xlim([-6.5 6.5])
        ylim([-1.5 2.6])
        xlabel('x (cm)', 'Interpreter', 'latex', 'FontWeight','bold')
        ylabel('y (cm)', 'Interpreter', 'latex', 'FontWeight','bold')
        lgd = legend(['Energy layer $E_{nominal}=' num2str(e(cnt)*1000) '$ MeV'], ...
            'Target area', 'Location', 'NorthWest', 'Interpreter', 'latex');
        title(lgd, '\bf{Scan spot map w/ Meterset weights}')
        text(ScanSpotPositionMap{cnt}(:,2) + 0.1, ...
            ScanSpotPositionMap{cnt}(:,1) + 0.1, ...
            cellstr(num2str(round(ScanSpotMetersetWeights{cnt},4))), ...
            'FontSize', 9);
        set(gca, 'FontSize', 12)
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
FWHM_x_mrad   = (2*sqrt(2*log(2))).*sigma_x_mrad;
FWHM_y_mrad   = (2*sqrt(2*log(2))).*sigma_y_mrad;

% Use 'real' instead of nominal beam energy
e = e_real;

% [e'.*1000 sigma_e.*1000 sigma_x.*10 sigma_y.*10]
% mean(sigma_e)*1000

%% Define energy interval for each proton beamlet
% Assumption: mean(diff(e)) <- not always true!!
% What are proton beamlet energy range in FLUKA!??!
% What is e_min and e_max actually?! i.e. how to find e_min and e_max from
% 83-107 MeV, n = 8
% 84.056-107.7720 Mev, sigma = 1.178-1.104 MeV, n = 8
%  -0.3025;
% e_sigma or e_FWHM?!
% e_min = [e(1)-mean(sigma_e) e(1:end-1)]; % mean(diff(e))
% e_max = e;

%% Write to file
% fid = fopen(fullfile(destpath, ...
%     ['SOBP_khi' num2str(khi*100) '_p' num2str(p*100) '.txt']), 'wt');
% fprintf(fid, '%.6f  %.6f   %.6f\n', [e_min; e_max; w]);
% fclose(fid);

fid = fopen(fullfile(destpath, ['SOBP_E0' num2str(e(end)*10^6) ...
    '_chi' num2str(chi*100) '_p' num2str(p*100) '.dat']), 'wt');
fprintf(fid, ['E[GeV]\t\tsigma_E[GeV]\tw_norm\t' ...
    'fwhm_x[cm]\tfwhm_y[cm]\tfwhm_x[mrad]\tfwhm_y[mrad]\n']);
fprintf(fid, ['%.6f\t%.6f\t %.6f\t %.6f\t %.6f\t %.6f\t%.6f\n'], ...
    [e; sigma_e; w; FWHM_x_cm; FWHM_y_cm; FWHM_x_mrad; FWHM_y_mrad]);
fclose(fid);

%% Analytical BPs

% z = linspace(0, 10, 1000);  % depth in cm
% alpha_prime_H20 = 8.7E-8;   % in GeV^2/cm
% sigma_mono  = sigma_mono(alpha_prime, alpha, p0, R0);
% sigma_tot   = sigma_tot(sigma_mono, sigma_E0, alpha, p0, E0);
%
% D_H20_hat   = Phi0
% D_H20       =

% D_SOBP = (358/(chi^0.435*E0^0.77))*1.602176462*10^(-10); % in Gy*cm2
D_SOBP          = 0.085; % 0.095; % in nGy/primary

%% FLUKA MC BPs

filelist = getAllFiles(sourcepath_SOBP{2});

for i = 1:length(filelist)

    [depth, dose] = getDepthDose(filelist{i});

    % Convert to Gy
    dose = dose .* 1.602176462*10^(-7) * 10^9; % in nGy
    % Normalize to max
    % Dose = (Dose ./ max(Dose(:))) * 100;

    % Compute SOBP
    if i == 1
        dose_SOBP = zeros(length(dose), 1);
    end
    dose_SOBP = dose_SOBP + w(i).*dose;

    % Legend text
    lgdstr{i} = ['$E=' num2str(round(e(i)*1000,1)) '$ MeV, ', ...
        '$\sigma_E=' num2str(round(sigma_e(i)*1000,1)) '$ MeV, ', ...
        '$\sigma_x=' num2str(round(sigma_x_cm(i)*10,2)) '$ mm, ', ...
        '$\sigma_y=' num2str(round(sigma_y_cm(i)*10,2)) '$ mm'];

    % Store beamlet values
    DepthDose{i}.depth          = depth;
    DepthDose{i}.dose           = dose;
    DepthDose{i}.w              = w(i);
    DepthDose{i}.dose_w         = w(i).*dose;
    DepthDose{i}.e              = [num2str(e(i)*1000) ' MeV'];
    DepthDose{i}.sigma_e        = [num2str(sigma_e(i)*1000) ' MeV'];
    DepthDose{i}.sigma_x_mm     = [num2str(sigma_x_cm(i)*10) ' mm'];
    DepthDose{i}.sigma_y_mm     = [num2str(sigma_y_cm(i)*10) ' mm'];

    % Plot
    figure(1)
    hold on
    plot(DepthDose{i}.depth, DepthDose{i}.dose, lineSpec{i}, ...
        'LineWidth', 1)
    hold off
    xlim([0 inf])
    ylim([0 inf])
    xlabel('Depth (cm)', 'Interpreter', 'latex', 'FontWeight','bold')
    ylabel('Dose (nGy per primary)', 'Interpreter', 'latex', 'FontWeight','bold') % nGy per primary, particles/cm^2 per primary
    legend(lgdstr, 'Location', 'NorthWest', 'Interpreter', 'latex');
    title(['FLUKA Monte Carlo (MC) simulated pristine Bragg peaks (BPs), ' ...
        'Jette D. \& Chen W., (2011)'], 'Interpreter', 'latex')
    set(gca, 'FontSize', 12)

end

lgdstr{end+1}   = 'SOBP; experimental DCPT beam model'; % 'SOBP; Jette D. \& Chen W., (2011)';
D_SOBP_const    = repmat(D_SOBP, 1, length(depth));

% Plot SOBP
figure(2)
hold on
for i = 1:length(DepthDose)

    plot(DepthDose{i}.depth, DepthDose{i}.dose_w, lineSpec{i}, ...
        'LineWidth', 1)

    % Update legend text
    lgdstr{i} = ['$E=' num2str(round(e(i)*1000,1)) '$ MeV, ', ...
        '$\sigma_E=' num2str(round(sigma_e(i)*1000,1)) '$ MeV, ', ...
        '$\sigma_x=' num2str(round(sigma_x_cm(i)*10,2)) '$ mm, ', ...
        '$\sigma_y=' num2str(round(sigma_y_cm(i)*10,2)) '$ mm, ', ...
        '$w=' num2str(round(w(i)*100,2)) '\%$'];

end
plot(depth, dose_SOBP, lineSpec{i+1}, 'LineWidth', 1)
% plot(depth, D_SOBP_const, lineSpec{i+2}, 'LineWidth', 1)
% lgdstr{i+2} = '$D_{SOBP}$; Jette D. \& Chen W., (2011)';
hold off
xlim([0 inf])
ylim([0 0.11])
xlabel('Depth (cm)', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Dose (nGy per primary)', 'Interpreter', 'latex', 'FontWeight', 'bold')
legend(lgdstr, 'Location', 'NorthWest', 'Interpreter', 'latex');
title(titlestr, 'Interpreter', 'latex')
set(gca, 'FontSize', 12)

%% Beamlet weight optimization

Dose_array = [DepthDose{1}.dose, DepthDose{2}.dose, DepthDose{3}.dose, ...
    DepthDose{4}.dose, DepthDose{5}.dose, DepthDose{6}.dose, ...
    DepthDose{7}.dose, DepthDose{8}.dose];
[m, n] = size(Dose_array);

BP = 5.5 <= depth' & depth' <= 8.5;

% (BP,:)   (BP)
w_opt = lsqlin(Dose_array, D_SOBP_const, [], [], ones(1,n), ...
    1, zeros(n,1), ones(n,1));

% Plot SOBP
figure(3)
hold on
for i = 1:length(DepthDose)

    % Compute SOBP
    if i == 1
        dose_SOBP_opt = zeros(length(m), 1);
    end
    dose_SOBP_opt = dose_SOBP_opt + w_opt(i).*DepthDose{i}.dose;

    % Store beamlet values
    DepthDose{i}.w_opt              = w_opt(i);
    DepthDose{i}.dose_w_opt         = w_opt(i).*DepthDose{i}.dose;

    % Update legend text
    lgdstr{i} = ['$E=$' num2str(round(e(i)*1000,1)) ' MeV, ', ...
        '$\sigma_E=$' num2str(round(sigma_e(i)*1000,1)) ' MeV, ', ...
        '$\sigma_x=$' num2str(round(sigma_x_cm(i)*10,2)) ' mm, ', ...
        '$\sigma_y=$' num2str(round(sigma_y_cm(i)*10,2)) ' mm, ', ...
        '$w_{opt}=' num2str(round(w_opt(i)*100,2)) '\%$'];

    % Plot
    plot(DepthDose{i}.depth, DepthDose{i}.dose_w_opt, lineSpec{i}, ...
        'LineWidth', 1)

end
plot(depth, dose_SOBP_opt, lineSpec{i+1}, 'LineWidth', 1)
hold off
lgdstr{length(DepthDose)+1} = 'Optimized SOBP by CLLS';
xlim([0 inf])
ylim([0 0.11])
xlabel('Depth (cm)', 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('Dose (nGy per primary)', 'Interpreter', 'latex', 'FontWeight', 'bold')
legend(lgdstr, 'Location', 'NorthWest', 'Interpreter', 'latex');
title(titlestr, 'Interpreter', 'latex')
set(gca, 'FontSize', 12)


%%
