clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path_calib      = 'C:\Users\delmo\Desktop\EBT3 - 14122022 - PROCCA\Calibration';
path_ctrl       = 'C:\Users\delmo\Desktop\EBT3 - 14122022 - PROCCA\Control';
path_bckg       = 'C:\Users\delmo\Desktop\EBT3 - 14122022 - PROCCA\Background';
path_Open       = 'C:\Users\delmo\Desktop\EBT3 - 14122022 - PROCCA\Dose Measurement\Open';
path_Collimated = 'C:\Users\delmo\Desktop\EBT3 - 14122022 - PROCCA\Dose Measurement\Collimated';
path_Control    = 'C:\Users\delmo\Desktop\EBT3 - 14122022 - PROCCA\Dose Measurement\Control';

path_dest = 'C:\Users\delmo\Desktop\GRID\Results_v6';

%% Variables

% Scanning resolution (in dpi), bit depth and ROI size (in mm)
dpi             = 300; % [1200 300 1200 300];
px_size         = (2.54/dpi);     % in cm/pixel
bit             = 48;
bitperchannel   = bit/3;
ROI_size_mm     = [4 4];

% EBT3 film dose (in Gy)
dose = repmat([0 0.5 10 1 2 5]', 1, 6)';
dose = dose(:)';

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = [1 460]; % [1  570];
y_range_px = [12 590]; % ;
% x_range_px = [1 472]; % [1  570];
% y_range_px = [1 354]; %

% Plot spesifications
channel             = {'red', 'green', 'blue', 'gray'};
markerspec          = {'ro', 'gx', 'bs', 'kd'};
linespec            = {'r-', 'g-', 'b-', 'k-'};
lgdstr_profiles     = {'Red channel', 'Green channel', 'Blue channel', ...
    'Grayscale image'};

%% Read EBT3 films

EBT3_calib              = getEBT3struct(path_calib, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_ctrl               = getEBT3struct(path_ctrl, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
img_bckg                = getEBT3struct(path_bckg, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_Open               = getEBT3struct(path_Open, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_Collimated         = getEBT3struct(path_Collimated, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_Control   = getEBT3struct(path_Control, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);

%% Loop over controls (0 Gy) films and opaque (background) images

for j = 1:length(channel)

    % Initialization per channel
    PVvec_bckg              = [];
    sigmavec_PV_bckg        = [];
    PVvec_ctrl              = [];
    sigmavec_PV_ctrl        = [];
    netODvec_ctrl           = [];
    sigmavec_netOD_ctrl     = [];

    % Loop over background images and vectorize their PV values
    for i = 1:length(img_bckg)
        PVvec_bckg         = [PVvec_bckg img_bckg{i}.(channel{j}).PV];
        sigmavec_PV_bckg   = [sigmavec_PV_bckg img_bckg{i}.(channel{j}).sigma];
    end

    % Estimate mean PV and SD for background images
    [bckg.(channel{j}).PV_avg, bckg.(channel{j}).sigma_PV_avg] = ...
        estimateAverageI(PVvec_bckg, sigmavec_PV_bckg);

    % Loop over control images and vectorize their netOD
    for i = 1:length(EBT3_ctrl)

        %         netOD_ctrl = EBT3_ctrl{i}.(channel{j}).OD - EBT3_ctrl{i}.(channel{j}).OD;
        netOD_ctrl = calculateNetOD(EBT3_ctrl{i}.(channel{j}).PV, ...
            EBT3_ctrl{i}.(channel{j}).PV, bckg.(channel{j}).PV_avg);
        sigma_netOD_ctrl = calculateSigmaNetOD( ...
            EBT3_ctrl{i}.(channel{j}).PV, EBT3_ctrl{i}.(channel{j}).PV, ...
            bckg.(channel{j}).PV_avg, EBT3_ctrl{i}.(channel{j}).sigma, ...
            EBT3_ctrl{i}.(channel{j}).sigma, bckg.(channel{j}).sigma_PV_avg);

        PVvec_ctrl              = [PVvec_ctrl EBT3_ctrl{i}.(channel{j}).PV];
        sigmavec_PV_ctrl        = [sigmavec_PV_ctrl EBT3_ctrl{i}.(channel{j}).sigma];
        netODvec_ctrl           = [netODvec_ctrl netOD_ctrl];
        sigmavec_netOD_ctrl     = [sigmavec_netOD_ctrl sigma_netOD_ctrl];

    end

    % Estimate and store mean PV and SD for control images
    [ctrl.(channel{j}).PV_avg, ctrl.(channel{j}).sigma_PV_avg] = ...
        estimateAverageI(PVvec_ctrl, sigmavec_PV_ctrl);

    % Estimate and store mean netOD and SD for control images
    [ctrl.(channel{j}).netOD_avg, ctrl.(channel{j}).sigma_netOD_avg] = ...
        estimateAverageI(netODvec_ctrl, sigmavec_netOD_ctrl);

    % Store netOD and SD for control images
    ctrl.(channel{j}).netOD         = netODvec_ctrl;
    ctrl.(channel{j}).sigma_netOD   = sigmavec_netOD_ctrl;

end

%% Loop over calibration (> 0 Gy) films

for j = 1:length(channel)

    %     ODvec         = [];
    netODvec        = ctrl.(channel{j}).netOD; % [];
    sigmavec_netOD  = ctrl.(channel{j}).sigma_netOD; % [];

    for i = 1:length(EBT3_calib)
        %         ODvec = [ODvec EBT3_calib{i}.(channel{j}).OD];
        netOD = calculateNetOD(ctrl.(channel{j}).PV_avg, ...
            EBT3_calib{i}.(channel{j}).PV, bckg.(channel{j}).PV_avg);
        sigma_netOD = calculateSigmaNetOD(ctrl.(channel{j}).PV_avg, ...
            EBT3_calib{i}.(channel{j}).PV, bckg.(channel{j}).PV_avg, ...
            ctrl.(channel{j}).sigma_PV_avg, ...
            EBT3_calib{i}.(channel{j}).sigma, ...
            bckg.(channel{j}).sigma_PV_avg);

        netODvec = [netODvec netOD];
        sigmavec_netOD = [sigmavec_netOD sigma_netOD];
    end

    calib.(channel{j}).netOD        = netODvec;
    calib.(channel{j}).sigma_netOD  = sigmavec_netOD;

end

%% Functional form regression fit

for j = 1:length(channel)

    RsquaredAdjusted_opt = 0;

    % Search space for the model parameter n, while keeping it fixed as a
    % constant
    for n = 0:0.1:5

        % Design matrix
        X = [calib.(channel{j}).netOD' (calib.(channel{j}).netOD').^n];

        % Fit linear model
        mdl = fitlm(X, dose', 'Intercept', false); % (N_control+1)

        % Store fitted model with best adjusted R-squared
        if mdl.Rsquared.Adjusted > RsquaredAdjusted_opt
            RsquaredAdjusted_opt = mdl.Rsquared.Adjusted;
            fit.(channel{j}).mdl = mdl;
            fit.(channel{j}).n   = n;
            fit.(channel{j}).RsquaredAdjusted = mdl.Rsquared.Adjusted;
            calib.(channel{j}).mdl = mdl;
            calib.(channel{j}).n   = n;
        end

    end

end

% Plot
%     figure();
hold on
for j = 1:length(channel)

    % Interpolate model fit to a finer grid
    netOD_interp = linspace(min(calib.(channel{j}).netOD), ...
        max(calib.(channel{j}).netOD), 10000);
    dose_interp = model1(netOD_interp, ...
        fit.(channel{j}).mdl.Coefficients.Estimate(1), ...
        fit.(channel{j}).mdl.Coefficients.Estimate(2), ...
        fit.(channel{j}).n);

    % Plot
    h_data(j) = plot(dose, calib.(channel{j}).netOD, markerspec{j});
    h_fit(j)  = plot(dose_interp, netOD_interp, linespec{j});

end
%     xlabel('Dose (Gy)')
%     ylabel('\it{netOD}') % 'Transmittance \it{T}' % '\it{OD}'
%     xlim([min(dose)-0.5 max(dose)+0.5])
%     legend([h_data(1), h_data(2), h_data(3), h_data(4), h_fit(4), ...
%         h_data(5), h_data(6), h_data(7), h_data(8), h_fit(8)], ...
%         'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
%         'Data, grayscale image', 'Fit (\it{a\cdotnetOD+b\cdotnetOD^n})', ...
%         'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
%         'Data, grayscale image', 'Fit (\it{a\cdotnetOD+b\cdotnetOD^n})', ...
%         'Location', 'NorthWest')
%     grid on
%     set(gca, 'FontSize', 16)
%     hold off

xlabel('Dose (Gy)')
ylabel('\it{netOD}') % 'Transmittance \it{T}' % '\it{OD}'
xlim([min(dose)-0.5 max(dose)+0.5])
legend([h_data(1), h_data(2), h_data(3), h_data(4), h_fit(4)], ...
    'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
    'Data, grayscale image', 'Fit (\it{a\cdotnetOD+b\cdotnetOD^{n}})', ...
    'Location', 'NorthWest')
grid on
set(gca, 'FontSize', 16)
hold off

%% Convert experimental EBT3 films into dose by using the calibration

% clc
% close all

channel_opt     = 'red';
pattern         = {'Open', 'Collimated', 'Control'};

center_EBT3     = 290;
x1_EBT3         = center_EBT3 - 10;
x2_EBT3         = center_EBT3 + 10;
significance    = 0.0001; % 0.000001 is too unstable for CI bands

h_profilevec = [];
dose_EBT3    = {};
profile      = {};

for i = 1:length(pattern) % 2

    img_dose = [];
    for j = 1:length(eval(sprintf('EBT3_%s', pattern{i})))

        temp_img_dose = convertEBT3netODtoDose( ...
            eval(sprintf('EBT3_%s', pattern{i})), ...
            channel_opt, fit, ctrl, bckg, 0);
        img_dose = cat(3, img_dose, temp_img_dose);

    end
    struct.(pattern{i}).img_dose = img_dose;

    %     temp_img_dose = convertEBT3netODtoDose( ...
    %         eval(sprintf('EBT3_%s', pattern{i})), ...
    %         channel_opt, fit, ctrl, bckg, 0);
    %     img_dose = cat(3, img_dose, temp_img_dose);
    %     struct.(pattern{i}).img_dose = img_dose;

    doseprofiles = [];
    if strcmp('Open', pattern{i})
        dose_EBT3.(pattern{i}).dosevalues        = [];
    else
        dose_EBT3.(pattern{i}).valley.dosevalues = [];
        dose_EBT3.(pattern{i}).peak.dosevalues   = [];
    end

    figure(100)
    hold on
    for j = 1:size(struct.(pattern{i}).img_dose, 3)

        % Dose profiles
        temp_img = struct.(pattern{i}).img_dose(1:end, x1_EBT3:x2_EBT3, j);
        temp_img(temp_img <= 0) = 0;

        temp_doseprofile    = mean(temp_img, 2);
        position            = (1:length(temp_doseprofile)) .* px_size;

        % Store dose profile
        doseprofiles = [doseprofiles; temp_doseprofile'];

        % Plot 1D profile across estimated EBT3 dose maps
        h = plot(position, temp_doseprofile, linespec{i}, 'LineWidth', 1.0);
        %     yline(dose_valley, '--', 'Valley')
        %     yline(dose_peak, '--', 'Peak')
        %     yline(dose_open, '-.', 'Open')

    end
    h_profilevec = [h_profilevec h];
    lgdstr_profiles{i} = [strrep(pattern{i}, '_', ' ') ', ($N=6$ films)'];
    xlabel('Position (cm)')
    ylabel('Dose (Gy)')
    xlim([0.375 4.2])
    ylim([0 10])
    legend(h_profilevec, lgdstr_profiles, 'Location', 'NorthEast', ...
        'Interpreter', 'LaTeX')
    grid on
    set(gca, 'FontSize', 20)
    hold off

    if strcmp('Open', pattern{i})
        [dose_EBT3.(pattern{i}).doseavg, dose_EBT3.(pattern{i}).CI] = ...
            estimateCI(dose_EBT3.(pattern{i}).dosevalues, significance);
    else
        [dose_EBT3.(pattern{i}).valley.doseavg, ...
            dose_EBT3.(pattern{i}).valley.CI] = ...
            estimateCI(dose_EBT3.(pattern{i}).valley.dosevalues, significance);
        [dose_EBT3.(pattern{i}).peak.doseavg, ...
            dose_EBT3.(pattern{i}).peak.CI] = ...
            estimateCI(dose_EBT3.(pattern{i}).peak.dosevalues, significance);
    end

    % Estimate 95% CI bands for the dose profiles
    profile.(pattern{i}).Dose = estimateCIProfileBand( ...
        doseprofiles, position, significance);

    % Normalization
    if i == 1
        norm = mean(profile.(pattern{i}).Dose.avgprofile);
    end
    profile.(pattern{i}).Dose.avgprofile = ...
        profile.(pattern{i}).Dose.avgprofile ./ norm;
    profile.(pattern{i}).Dose.CI95_value = ...
        profile.(pattern{i}).Dose.CI95_value ./ norm;

end

%% Plot EBT3 2D dose maps

for i = 1:length(pattern)

    plot2Dmap(struct.(pattern{i}).img_dose, px_size, [0 8], 'Dose (Gy)')

end
