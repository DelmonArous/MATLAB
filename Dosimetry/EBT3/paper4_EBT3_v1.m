clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path_calib          = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Calibration_v1';
path_ctrl           = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Control_v1';
path_bckg           = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Background';
path_Open           = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Dose Measurement\Open_v1';
path_GRID_Stripes   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Dose Measurement\Stripes_v1';
path_GRID_Dots      = 'C:\Users\delmo\Desktop\Jacob\EBT3\131021 - Holes\Dose Measurements\Holes_v1';
% path_assay        = 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation_v6';
path_assay          = 'C:\Users\delmo\Desktop\Jacob\Segmentation Binary Masks_v6';
path_bwflask        = 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks_v6';

path_dest = 'C:\Users\delmo\Desktop\GRID\Results_v6';

%% Variables

% Scanning resolution (in dpi), bit depth and ROI size (in mm)
dpi             = 300; % [1200 300 1200 300];
px_size         = (2.54/dpi);     % in cm/pixel
bit             = 48;
bitperchannel   = bit/3;
ROI_size_mm     = [4 4];

% EBT3 film dose (in Gy)
dose = repmat([0 0.1 0.2 0.5 10 1 2 5]', 1, 8)';
dose = dose(:)';

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = [1 470]; % [1  570];
y_range_px = [10 500]; % ;
% x_range_px = [1 472]; % [1  570];
% y_range_px = [1 354]; %

% Plot spesifications
channel             = {'red', 'green', 'blue', 'gray'};
% markerspec          = {'ro', 'gx', 'bs', 'kd'};
markerspec          = {{'ro', 'go', 'bo', 'ko'}, {'rs', 'gs', 'bs', 'ks'}};
% linespec            = {'r-', 'g-', 'b-', 'k-'};
linespec            = {{'r-', 'g-', 'b-', 'k-'}, {'r--', 'g--', 'b--', 'k--'}};
lgdstr_profiles     = {'Red channel', 'Green channel', 'Blue channel', ...
    'Grayscale image'};

%% Read EBT3 films

EBT3_calib          = getEBT3struct(path_calib, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_ctrl           = getEBT3struct(path_ctrl, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
img_bckg            = getEBT3struct(path_bckg, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_Open           = getEBT3struct(path_Open, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_GRID_Stripes   = getEBT3struct(path_GRID_Stripes, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);
EBT3_GRID_Dots      = getEBT3struct(path_GRID_Dots, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 1);

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

% Differential response
ind_calib   = {[1:10 12:13 15:18 23:36 41:46 51:61 63:64], ...
    [1:8 11 14 19:22 37:40 47:50 62]};
response    = {'high', 'low'};

for i = 1:length(ind_calib)

    for j = 1:length(channel)

        RsquaredAdjusted_opt = 0;

        % Search space for the model parameter n, while keeping it fixed as a
        % constant
        for n = 0:0.1:5

            % Design matrix
            X = [calib.(channel{j}).netOD(ind_calib{i})' ...
                (calib.(channel{j}).netOD(ind_calib{i})').^n];

            % Fit linear model
            mdl = fitlm(X, dose(ind_calib{i})', 'Intercept', false); % (N_control+1)

            % Store fitted model with best adjusted R-squared
            if mdl.Rsquared.Adjusted > RsquaredAdjusted_opt
                RsquaredAdjusted_opt = mdl.Rsquared.Adjusted;
                fit.(response{i}).(channel{j}).mdl = mdl;
                fit.(response{i}).(channel{j}).n = n;
                fit.(response{i}).(channel{j}).RsquaredAdjusted = mdl.Rsquared.Adjusted;
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
        netOD_interp = linspace( ...
            min(calib.(channel{j}).netOD(ind_calib{i})), ...
            max(calib.(channel{j}).netOD(ind_calib{i})), 10000);
        dose_interp = model1(netOD_interp, ...
            fit.(response{i}).(channel{j}).mdl.Coefficients.Estimate(1), ...
            fit.(response{i}).(channel{j}).mdl.Coefficients.Estimate(2), ...
            fit.(response{i}).(channel{j}).n);

        % Plot
        h_data((i-1)*length(channel) + j)   = plot(dose(ind_calib{i}), ...
            calib.(channel{j}).netOD(ind_calib{i}), markerspec{i}{j});
        h_fit((i-1)*length(channel) + j)    = plot(dose_interp, ...
            netOD_interp, linespec{i}{j});

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

end
xlabel('Dose (Gy)')
ylabel('\it{netOD}') % 'Transmittance \it{T}' % '\it{OD}'
xlim([min(dose)-0.5 max(dose)+0.5])
legend([h_data(1), h_data(2), h_data(3), h_data(4), h_fit(4), ...
    h_data(5), h_data(6), h_data(7), h_data(8), h_fit(8)], ...
    'Data (high), red channel', 'Data (high), green channel', ...
    'Data (high), blue channel', 'Data (high), grayscale image', ...
    'Fit (high) (\it{a_h\cdotnetOD+b_h\cdotnetOD^{n_h}})', ...
    'Data (low), red channel', 'Data (low), green channel', ...
    'Data (low), blue channel', 'Data (low), grayscale image', ...
    'Fit (low) (\it{a_l\cdotnetOD+b_l\cdotnetOD^{n_l}})', 'Location', 'NorthWest')
% legend([h_data(1), h_data(2), h_data(3), h_data(4), h_fit(4)], ...
%     'Data, red channel', 'Data, green channel', ...
%     'Data, blue channel', 'Data, grayscale image', ...
%     'Fit (\it{a\cdotnetOD+b\cdotnetOD^{n}})', 'Location', 'NorthWest')
grid on
set(gca, 'FontSize', 16)
hold off

%% Convert experimental EBT3 films into dose by using the calibration

% clc
close all

channel_opt = 'red';

pattern         = {'Open', 'GRID_Stripes', 'GRID_Dots'};
linespec        = {'-r', '-g', '-b'};

center_EBT3     = 282;
x1_EBT3         = center_EBT3 - 10;
x2_EBT3         = center_EBT3 + 10;
significance    = 0.0001; % 0.000001 is too unstable for CI bands

% Valley and peak positions
pos                     = {};
pos.GRID_Stripes.valley = [60:160 250:350 440:540];
pos.GRID_Stripes.peak   = [185:230 375:410 560:600];
pos.GRID_Dots.valley    = [45:55 140:360 450:675];
pos.GRID_Dots.peak      = [80:110 390:420];
pos.Open                = 40:700; % pos.GRID_dots.valley(1):pos.GRID_dots.peak(end);

h_profilevec = [];
dose_EBT3    = {};
profile      = {};

% Differential response (high, low)
ind.Open         = {[1 4:5 10:14],  [2:3 6:9 15:16]};
ind.GRID_Stripes = {[1:7 12 14:16], [8:11 13]};
ind.GRID_Dots    = {[7:10 12:15], [1:6 11 16]};

for i = 1:length(pattern) % 2

    img_dose = [];
    for j = 1:length(ind.(pattern{i}))

        temp_img_dose = convertEBT3netODtoDose( ...
            eval(sprintf('EBT3_%s', pattern{i})), ...
            channel_opt, fit.(response{j}), ctrl, bckg, ind.(pattern{i}){j});
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

        % Save peak and valley dose for each profile of a pattern
        if strcmp('Open', pattern{i})
            dose_EBT3.(pattern{i}).dosevalues = [ ...
                dose_EBT3.(pattern{i}).dosevalues ...
                mean(temp_doseprofile(pos.(pattern{i})))];
        else
            dose_EBT3.(pattern{i}).valley.dosevalues = [ ...
                dose_EBT3.(pattern{i}).valley.dosevalues ...
                mean(temp_doseprofile(pos.(pattern{i}).valley))];
            dose_EBT3.(pattern{i}).peak.dosevalues = [ ...
                dose_EBT3.(pattern{i}).peak.dosevalues ...
                mean(temp_doseprofile(pos.(pattern{i}).peak))];
        end

        % Store dose profile
        doseprofiles = [doseprofiles; temp_doseprofile'];

        % Plot 1D profile across estimated EBT3 dose maps
        h = plot(position, temp_doseprofile, linespec{i}, 'LineWidth', 1.0);
        %     yline(dose_valley, '--', 'Valley')
        %     yline(dose_peak, '--', 'Peak')
        %     yline(dose_open, '-.', 'Open')

    end
    h_profilevec = [h_profilevec h];
    lgdstr_profiles{i} = [strrep(pattern{i}, '_', ' ') ' ($N=16$ films)'];
    xlabel('Position (cm)')
    ylabel('Dose (Gy)')
    xlim([0.375 5.9])
    ylim([0 7])
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
    profile.(pattern{i}).Dose = estimateCIBand(doseprofiles, position, significance);

    % Normalization
    if i == 1
        norm = mean(profile.(pattern{i}).Dose.y_avg);
    end
    profile.(pattern{i}).Dose.y_avg = ...
        profile.(pattern{i}).Dose.y_avg ./ norm;
    profile.(pattern{i}).Dose.CI95_value = ...
        profile.(pattern{i}).Dose.CI95_value ./ norm;

end

%% 95% CI for the dose profiles

lgdstr = {'EBT3 mean profile, Open', 'EBT3 95\% CI, Open', ...
    'EBT3 mean profile, Stripes', 'EBT3 95\% CI, Stripes', ...
    'EBT3 mean profile, Dots', 'EBT3 95\% CI, Dots'};
color_CI95     = {'red', 'green', 'blue'};
facecolor_CI95 = {[1.0 0.8 0.8], [0.78 1.0 0.46], [0.3010 0.7450 0.9330]};

h_vec = [];

for i = 1:length(pattern)

    figure(200)
    hold on
    h_95CIband = fill(profile.(pattern{i}).Dose.CI95_x, ...
        profile.(pattern{i}).Dose.CI95_value, color_CI95{i}, ...
        'FaceColor', facecolor_CI95{i}, 'EdgeColor', 'none');
    h_95CImean = plot(profile.(pattern{i}).Dose.x, ...
        profile.(pattern{i}).Dose.y_avg, ...
        linespec{i}, 'LineWidth', 1.0);
    hold off
    h_vec = [h_vec h_95CImean h_95CIband];

end
legend(h_vec, lgdstr, 'Location', 'NorthEast', 'NumColumns', 3, ...
    'Interpreter', 'LaTeX')
xlabel('Position (cm)')
ylabel('Normalized Dose')
xlim([0.375 5.8])
ylim([0 1.2])
grid on
set(gca, 'FontSize', 22)

%% Compute peak distance image

% Peak distance coordinates for computation of peak distance image
peak_xycoord.GRID_Dots = [282 97; 127 251; 438 250; 283 404; ...
    130 562; 441 560; 286 714];
peak_xycoord.GRID_Stripes = [...
    (25:480)' repmat(39, 1, length(25:480))'; ...
    (25:480)' repmat(206, 1, length(25:480))'; ...
    (25:480)' repmat(392, 1, length(25:480))'; ...
    (55:417)' repmat(580, 1, length((55:417)))'];

% Peak area variable
peakarea.Control      = 0.0;
peakarea.GRID_Dots    = 0.05;
peakarea.GRID_Stripes = 0.33;
peakarea.Open         = 1.0;

for j = 1:length(pattern)

    % Store dose map
    temp_img_dose = sum(struct.(pattern{j}).img_dose, 3) ./ ...
        size(struct.(pattern{j}).img_dose, 3);
    temp_img_dose(temp_img_dose <= 0) = 0;
    temp_img_dose(isnan(temp_img_dose)) = 0;
    temp_img_dose(imag(temp_img_dose) ~= 0) = 0;

    % Binarize dose map
    temp_img = (temp_img_dose-min(temp_img_dose(:))) ./ ...
        abs(max(temp_img_dose(:))-min(temp_img_dose(:)));
    temp_bw_dose = imbinarize(temp_img, 'adaptive');
    temp_bw_dose = imdilate(temp_bw_dose, strel('disk', 5));
    temp_bw_dose = bwmorph(temp_bw_dose, 'thicken');
    temp_bw_dose = bwmorph(temp_bw_dose, 'majority');
    temp_bw_dose = bwmorph(temp_bw_dose, 'bridge');
    temp_bw_dose = imclose(temp_bw_dose, strel('disk', 60));
    temp_bw_dose = imfill(temp_bw_dose, 'holes');
    [y, x]       = find(temp_bw_dose == 1);

    % Compute peak distance image
    if any(strcmp(pattern{j}, {'GRID_Dots', 'GRID_Stripes'}))
        [img_peakdist, ~] = EstimatePeakDistanceMap(temp_bw_dose, ...
            [x y], peak_xycoord.(pattern{j}));
    else
        img_peakdist = zeros(size(temp_img_dose,1), size(temp_img_dose,2));
    end
    img_peakdist = px_size .* img_peakdist;

    % Compute peak area image
%     img_peakarea = peakarea.(pattern{j}) .* ones(size(temp_img_dose,1), ...
%         size(temp_img_dose,2));
%     img_peakarea = temp_bw_dose .* img_peakarea;

    % Store peak distance and area image
    struct.(pattern{j}).img_peakdist = img_peakdist;
%     struct.(pattern{j}).img_grad     = img_grad;
%     struct.(pattern{j}).img_peakarea = img_peakarea;

    % Plot EBT3 2D dose and peak distance map
    plot2Dmap(struct.(pattern{j}).img_dose, px_size, [0 6], ...
        'Dose (Gy)')
    plot2Dmap(struct.(pattern{j}).img_peakdist, px_size, ...
        [0 max(max(struct.(pattern{j}).img_peakdist))], ...
        'Peak distance (cm)')
%     plot2Dmap(struct.(pattern{j}).img_grad, px_size, ...
%         [0 max(max(struct.(pattern{j}).img_grad))], ...
%         'Dose gradient (cm/Gy)')
%     plot2Dmap(img_peakarea, px_size, [0 max(max(img_peakarea))], ...
%         'Peak area ratio')

end

%%
close all
clc

px_size_assay   = 2.54/1200;    % in cm/px, 1200 dpi
cells_seeded    = 10000;        % no. of seeded cells (30000 originally)

lgdstr = {'Mean predicted survival', '95\% CI predicted', ...
    'Mean observed survival', '95\% CI observed', 'RPD'};
y_limits_SF     = {[0.5 1.25],     [0.15 1.2],  [0.0 1.1]};
y_limits_RPE    = {[-40 40],       [-80 90],    [-170 250]};

% Band/quadrat size
dx_cm       = 0.10;                           % band size in cm
dx_px       = round(dx_cm/px_size_assay);     % band size in px
dxdy_cm     = [0.10 0.10];                    % quadrat size in cm
dxdy_px     = round(dxdy_cm/px_size_assay);   % quadrat size in px

model           = {'LQ', 'MLQ'};
GRID_region     = {'peak', 'valley'};
pattern         = {'Control', 'Dotted', 'Striped', 'Open'};
assay_pattern   = {'Control', 'GRID_Dots', 'GRID_Stripes', 'Open'};
assay_dose      = {'Dose2Gy', 'Dose5Gy', 'Dose10Gy'};
dose_str        = {'2 Gy', '5 Gy', '10 Gy'};
dose_assay      = repmat([2 5 10]', 1, 4)';
dose_assay      = dose_assay(:)';
[~, ind_sort]   = sort(dose_assay, 'ascend');
dose_scale      = [2/5 1 2];

% Peak and valley doses
for i = 1:length(assay_dose)

    temp_dose = dose_scale(i) * 5;

    for j = 2:length(assay_pattern)

        if strcmp(assay_pattern{j}, 'Open')
            irrad_dose.(assay_dose{i}).(assay_pattern{j}).dose = temp_dose;
            irrad_dose.(assay_dose{i}).(assay_pattern{j}).dose_CI95 = temp_dose .* [0.96 1.02];
        elseif strcmp(assay_pattern{j}, 'GRID_Dots')
            temp_region_dose      = temp_dose .* [0.70 0.10];
            temp_region_dose_95CI = temp_dose .* [0.66 0.74; 0.091 0.13];
            for k = 1:length(GRID_region)
                irrad_dose.(assay_dose{i}).(assay_pattern{j}).(GRID_region{k}).dose = temp_region_dose(k);
                irrad_dose.(assay_dose{i}).(assay_pattern{j}).(GRID_region{k}).dose_CI95 = temp_region_dose_95CI(k,:);
            end
        else
            temp_region_dose        = temp_dose .* [0.82 0.18];
            temp_region_dose_95CI   = temp_dose .* [0.78 0.85; 0.16 0.19];
            for k = 1:length(GRID_region)
                irrad_dose.(assay_dose{i}).(assay_pattern{j}).(GRID_region{k}).dose = temp_region_dose(k);
                irrad_dose.(assay_dose{i}).(assay_pattern{j}).(GRID_region{k}).dose_CI95 = temp_region_dose_95CI(k,:);
            end
        end

    end
end

% Peak distance coordinates for computation of peak distance variable
peak_xycoord.GRID_Dots = [1083 414; 468 1037; 1710 1024; 1090 1650; ...
    474 2282; 1720 2270; 1096 2890];
peak_xycoord.GRID_Stripes = [(145:2025)' repmat(238, 1, length(145:2025))'; ...
    (145:2025)' repmat(905, 1, length(145:2025))'; ...
    (145:2005)' repmat(1655, 1, length(145:2005))'; ...
    (335:1795)' repmat(2400, 1, length(335:1795))'];

% Translation vectors for manual image registration
translation_vec.GRID_Dots = {[-120 55], [5 78], [-30 73], [-22 -3], ...
    [-70 45], [10 48], [85 8], [-16 1], ...
    [-50 43], [0 69], [89 -5], [0 0]};
translation_vec.GRID_Stripes = {[0 67], [0 60], [0 -15], [0 55], ...
    [0 37], [0 44], [0 1], [0 -11], ...
    [0 51], [0 48], [0 16], [0 0]};
translation_dose_vec.GRID_Dots = [-41 37];
translation_dose_vec.GRID_Stripes = [96 85];

folderlist_assay    = getAllFolders(path_assay);
folderlist_bwflask  = getAllFolders(path_bwflask);

for i = 1:length(folderlist_assay)

    [~, pattern_name, ~] = fileparts(folderlist_assay{i});
    filelist_assay      = getAllFiles(folderlist_assay{i});
    filelist_bwflask    = getAllFiles(folderlist_bwflask{i});

    if i > 1
        temp_img_dose = sum(struct.(assay_pattern{i}).img_dose, 3) ./ ...
            size(struct.(assay_pattern{i}).img_dose, 3);
        temp_img_dose(temp_img_dose <= 0) = 0;
        temp_img_dose(isnan(temp_img_dose)) = 0;
        temp_img_dose(imag(temp_img_dose) ~= 0) = 0;
        temp_img_peakdist = struct.(assay_pattern{i}).img_peakdist;
    end

    for j = 1:length(filelist_assay)

        [~, filename, ~] = fileparts(filelist_assay{j});
        bw_assay_seg = logical(readmatrix(filelist_assay{j}));
        bw_flask     = logical(readmatrix(filelist_bwflask{j}));

        % Manual image registration
        if contains(pattern_name, ["Dots", "Stripes"]) % && j > 4
            bw_assay_seg = imtranslate(bw_assay_seg, translation_vec.(assay_pattern{i}){j});
            bw_flask     = imtranslate(bw_flask, translation_vec.(assay_pattern{i}){j});
        end

        stats       = regionprops(bw_assay_seg, 'Centroid');
        centroids   = round(cat(1, stats.Centroid));
        x_centroids = centroids(:,2);
        y_centroids = centroids(:,1);
        N_colonies  = length(stats);
        assay.(assay_pattern{i}).img{j}.filename    = filename;
        assay.(assay_pattern{i}).img{j}.bw_seg      = bw_assay_seg;
        assay.(assay_pattern{i}).img{j}.bw_flask    = bw_flask;
        assay.(assay_pattern{i}).img{j}.x_centroids = x_centroids;
        assay.(assay_pattern{i}).img{j}.y_centroids = y_centroids;
%         assay.(assay_pattern{i}).img{j}.N_colonies  = length(stats);
        if strcmp(assay_pattern{i}, 'Control')
            assay.(assay_pattern{i}).img{j}.PE          = ...
                (length(stats) * f(0)) / cells_seeded;
            assay.(assay_pattern{i}).img{j}.dose    = 0.0;
        else
            assay.(assay_pattern{i}).img{j}.dose    = dose_assay(j);
        end

        % 2D survival analysis
        n_quadrat_x = ceil(size(bw_flask,1)./dxdy_px(1));
        n_quadrat_y = ceil(size(bw_flask,2)./dxdy_px(2));
        dummy_zeros_x = zeros(n_quadrat_x*dxdy_px(1)-size(bw_flask, 1), ...
            size(bw_flask, 2));
        dummy_zeros_y = zeros(n_quadrat_x*dxdy_px(1), ...
            n_quadrat_y*dxdy_px(2)-size(bw_flask, 2));
        bw_flask      = cat(1, bw_flask, dummy_zeros_x);
        bw_flask      = logical(cat(2, bw_flask, dummy_zeros_y));
        bw_assay_seg  = cat(1, bw_assay_seg, dummy_zeros_x);
        bw_assay_seg  = logical(cat(2, bw_assay_seg, dummy_zeros_y));
        bw_centroids  = false(size(bw_assay_seg,1), size(bw_assay_seg,2));
        bw_centroids(sub2ind(size(bw_assay_seg), ...
            x_centroids, y_centroids)) = 1;

        % Dividing the segmentation and flask mask into
        % (n_quadrat_y x n_quadrat_x) quadrats
        bw_flask_quadrat       = mat2quadrat(bw_flask, dxdy_px);
        bw_assay_seg_quadrat   = mat2quadrat(bw_assay_seg, dxdy_px);
        bw_centroids_quadrat   = mat2quadrat(bw_centroids, dxdy_px);
        temp_bw_flask_quadrat  = false(n_quadrat_x, n_quadrat_y);
        temp_img_count_quadrat = zeros(n_quadrat_x, n_quadrat_y);

        for k = 1:n_quadrat_x
            for l = 1:n_quadrat_y
                [y_temp, x_temp] = find(bw_centroids_quadrat{k,l});
                centroid_quadrat{k,l}.bw_centroid   = bw_centroids_quadrat{k,l};
                centroid_quadrat{k,l}.x_centroid    = x_temp;
                centroid_quadrat{k,l}.y_centroid    = y_temp;
                centroid_quadrat{k,l}.colony_count  = length(x_temp);
                temp_img_count_quadrat(k,l)         = centroid_quadrat{k,l}.colony_count;
                if ~any(bw_flask_quadrat{k,l}(:) == 0)
                    temp_bw_flask_quadrat(k,l) = 1;
                end
            end
        end
        
        % Store quadrat results
        assay.(assay_pattern{i}).img{j}.bw_flask             = bw_flask;
        assay.(assay_pattern{i}).img{j}.bw_seg               = bw_assay_seg;
        assay.(assay_pattern{i}).img{j}.bw_centroids         = bw_centroids;
        assay.(assay_pattern{i}).img{j}.n_quadrat_x          = n_quadrat_x;
        assay.(assay_pattern{i}).img{j}.n_quadrat_y          = n_quadrat_y;
        assay.(assay_pattern{i}).img{j}.bw_flask_quadrat     = temp_bw_flask_quadrat;
        assay.(assay_pattern{i}).img{j}.bw_seg_quadrat       = bw_assay_seg_quadrat;

        % Scale dose map accordingly to nominal dose
        if i == 1
            temp_img_dose_scaled = 0.0 .* temp_img_dose;
        else
            temp_img_dose_scaled = (dose_assay(j)/5) .* temp_img_dose;
        end

        % Interpolate dose maps to match spatial resolution of assay images
        [X,Y] = meshgrid(1:size(temp_img_dose_scaled,2), ...
            1:size(temp_img_dose_scaled,1));
        [X2,Y2] = meshgrid( ...
            1:px_size_assay/px_size:size(temp_img_dose_scaled,2), ...
            1:px_size_assay/px_size:size(temp_img_dose_scaled,1));
        interp_img_dose = interp2(X, Y, temp_img_dose_scaled, ...
            X2, Y2, 'spline');
        dummy_zeros_x = zeros(size(bw_assay_seg,1)-size(interp_img_dose,1), ...
            size(interp_img_dose,2));
        dummy_zeros_y = zeros(size(bw_assay_seg,1), ...
            size(bw_assay_seg,2)-size(interp_img_dose,2));
        interp_img_dose = cat(1, interp_img_dose, dummy_zeros_x);
        interp_img_dose = cat(2, interp_img_dose, dummy_zeros_y);

        if i == 1

            dose_quadrat           = zeros(n_quadrat_x, n_quadrat_y);
            peakdist_quadrat       = zeros(n_quadrat_x, n_quadrat_y);
            grad_quadrat           = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_count_quadrat = round(f(0) .* temp_img_count_quadrat);

            PV_dose_quadrat     = dose_quadrat(temp_bw_flask_quadrat);
            PV_peakdist_quadrat = peakdist_quadrat(temp_bw_flask_quadrat); 
            PV_grad_quadrat     = grad_quadrat(temp_bw_flask_quadrat);
            PV_count_quadrat    = temp_img_count_quadrat(temp_bw_flask_quadrat);

            % Store quadrat results
            assay.(assay_pattern{i}).img{j}.dose_avg             = 0;
            assay.(assay_pattern{i}).img{j}.N_colonies           = f(0)*N_colonies;
            assay.(assay_pattern{i}).img{j}.N_quadrat            = nnz(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.bw_quadrat           = temp_bw_flask_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose_quadrat     = dose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_peakdist_quadrat = dose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_grad_quadrat     = grad_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose             = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_peakdist         = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_grad             = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.PV.dose_quadrat      = PV_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.peakdist_quadrat  = PV_peakdist_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.grad_quadrat      = PV_grad_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.count_quadrat     = PV_count_quadrat;

        else

            % Interpolate peak distance and peak area maps
            interp_img_peakdist = interp2(X, Y, temp_img_peakdist, ...
                X2, Y2, 'spline');
            interp_img_peakdist = cat(1, interp_img_peakdist, dummy_zeros_x);
            interp_img_peakdist = cat(2, interp_img_peakdist, dummy_zeros_y);

            % Binarize dose map
            temp_interp_img_dose = (interp_img_dose-min(interp_img_dose(:))) ./ ...
                abs(max(interp_img_dose(:))-min(interp_img_dose(:)));
            bw_dose = imbinarize(temp_interp_img_dose, 'adaptive');
            bw_dose = imdilate(bw_dose, strel('disk', 5));
            bw_dose = bwmorph(bw_dose, 'thicken');
            bw_dose = bwmorph(bw_dose, 'majority');
            bw_dose = bwmorph(bw_dose, 'bridge');
            bw_dose = imclose(bw_dose, strel('disk', 60));
            bw_dose = imfill(bw_dose, 'holes');

            if j == 1 && strcmp(assay_pattern{i}, 'Open')

                [optimizer, metric] = imregconfig('multimodal');
                optimizer.InitialRadius = 0.009;
                optimizer.Epsilon = 1.5e-4;
                optimizer.GrowthFactor = 1.01;
                optimizer.MaximumIterations = 300;
                tform = imregtform(double(bw_dose), double(bw_flask), ...
                    'translation', optimizer, metric); % similarity translation
                Rfixed = imref2d(size(bw_flask));
                bw_dose = imwarp(bw_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_dose = imwarp(interp_img_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_peakdist = imwarp(interp_img_peakdist, tform, ...
                    'OutputView', Rfixed);

            elseif j > 1 && strcmp(assay_pattern{i}, 'Open')

                bw_dose = imwarp(bw_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_dose = imwarp(interp_img_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_peakdist = imwarp(interp_img_peakdist, tform, ...
                    'OutputView', Rfixed);

            elseif any(strcmp(assay_pattern{i}, {'GRID_Stripes', 'GRID_Dots'}))

                bw_dose = imtranslate(bw_dose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_dose = imtranslate(interp_img_dose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_peakdist = imtranslate(interp_img_peakdist, ...
                    translation_dose_vec.(assay_pattern{i}));

            end
            bw_dose = imerode(bw_dose, strel('disk', 40));
            bw_dose = bwmorph(bw_dose, 'shrink');
            bw_dose = bwmorph(bw_dose, 'thin');
            bw_dose = logical(bw_dose);

            % Locate peak and valley pixels (regions)
            if strcmp(assay_pattern{i}, 'Open')
                temp_bw_dose_peak = logical( ...
                    interp_img_dose >= floor(dose_assay(j)*0.96*10)/10 ...   % dose_assay(j)*0.96
                    & interp_img_dose <= ceil(dose_assay(j)*1.02*10)/10);    % dose_assay(j)*1.02
                temp_bw_dose_valley = false(size(interp_img_dose,1), ...
                    size(interp_img_dose,2));
            elseif strcmp(assay_pattern{i}, 'GRID_Stripes')
                temp_bw_dose_peak = logical( ...
                    interp_img_dose >= floor(dose_assay(j)*0.78*10)/10 ...   % dose_assay(j)*0.78
                    & interp_img_dose <= ceil(dose_assay(j)*0.85*10)/10);    % dose_assay(j)*0.85
                temp_bw_dose_valley = logical( ...
                    interp_img_dose >= floor(dose_assay(j)*0.16*10)/10 ...   % dose_assay(j)*0.16
                    & interp_img_dose <= ceil(dose_assay(j)*0.19*10)/10);    % dose_assay(j)*0.19
            elseif strcmp(assay_pattern{i}, 'GRID_Dots')
                temp_bw_dose_peak = logical( ...
                    interp_img_dose >= floor(dose_assay(j)*0.66*10)/10 ...   % dose_assay(j)*0.66
                    & interp_img_dose <= ceil(dose_assay(j)*0.74*10)/10);    % dose_assay(j)*0.74
                temp_bw_dose_valley = logical( ...
                    interp_img_dose >= floor(dose_assay(j)*0.091*10)/10 ...  % dose_assay(j)*0.091
                    & interp_img_dose <= ceil(dose_assay(j)*0.13*10)/10);    % dose_assay(j)*0.13
            end

            % Dividing the dose map, and corresponding binary mask, into
            % (n_quadrat_y x n_quadrat_x) quadrats
            bw_dose_quadrat         = mat2quadrat(bw_dose, dxdy_px);
            bw_dose_peak_quadrat    = mat2quadrat(temp_bw_dose_peak, dxdy_px);
            bw_dose_valley_quadrat  = mat2quadrat(temp_bw_dose_valley, dxdy_px);
            img_dose_quadrat        = mat2quadrat(interp_img_dose, dxdy_px);
            img_peakdist_quadrat    = mat2quadrat(interp_img_peakdist, dxdy_px);

            temp_bw_dose_quadrat        = false(n_quadrat_x, n_quadrat_y);
            temp_bw_dose_peak_quadrat   = false(n_quadrat_x, n_quadrat_y);
            temp_bw_dose_valley_quadrat = false(n_quadrat_x, n_quadrat_y);
            temp_img_dose_quadrat       = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_peakdist_quadrat   = zeros(n_quadrat_x, n_quadrat_y);

            for k = 1:n_quadrat_x
                for l = 1:n_quadrat_y
                    if ~any(bw_dose_quadrat{k,l}(:) == 0)
                        temp_bw_dose_quadrat(k,l) = 1;
                    end
                    if ~any(bw_dose_peak_quadrat{k,l}(:) == 0)
                        temp_bw_dose_peak_quadrat(k,l) = 1;
                    end
                    if ~any(bw_dose_valley_quadrat{k,l}(:) == 0)
                        temp_bw_dose_valley_quadrat(k,l) = 1;
                    end
                    temp_img_dose_quadrat(k,l) = mean(img_dose_quadrat{k,l}(:));
                    temp_img_peakdist_quadrat(k,l) = mean(img_peakdist_quadrat{k,l}(:));
                end
            end            
            temp_img_count_quadrat = round(f(temp_img_dose_quadrat) ...
                .* temp_img_count_quadrat);

            % Compute gradient image (cm/Gy)
%             interp_img_grad =  interp_img_peakdist ./ interp_img_dose;
%             temp_img_grad_quadrat = temp_img_peakdist_quadrat ./ temp_img_dose_quadrat;
            [interp_img_grad, ~]        = imgradient(interp_img_dose);
            [temp_img_grad_quadrat, ~]  = imgradient(temp_img_dose_quadrat);

%             temp_img_grad_quadrat(temp_bw_dose_quadrat == 0)    = 0;
%             temp_img_grad_quadrat(temp_img_grad_quadrat <= 0)   = 0;
%             temp_img_grad_quadrat(isnan(temp_img_grad_quadrat)) = 0;
%             temp_img_grad_quadrat(isinf(temp_img_grad_quadrat)) = 0;

            temp_bw                = bw_flask & bw_dose;
            temp_bw_peak           = temp_bw_dose_peak & temp_bw;
            temp_bw_valley         = temp_bw_dose_valley & temp_bw;
            temp_bw_quadrat        = temp_bw_flask_quadrat & temp_bw_dose_quadrat;
            temp_bw_peak_quadrat   = temp_bw_dose_peak_quadrat & temp_bw_quadrat;
            temp_bw_valley_quadrat = temp_bw_dose_valley_quadrat & temp_bw_quadrat;
%            temp_bw_quadrat = temp_bw_dose_peak_quadrat | temp_bw_dose_valley_quadrat; % NEW!!!

%             min_img_peakdist_quadrat = prctile(temp_img_peakdist_quadrat(:), 25); % min(temp_img_peakdist_quadrat(:));
%             max_img_dose_quadrat     = prctile(temp_img_dose_quadrat(:), 75); % max(temp_img_dose_quadrat(:));
%             min_img_grad_quadrat     = prctile(temp_img_grad_quadrat(:), 25); % min(temp_img_grad_quadrat(:)); % R/D
%             max_img_grad_quadrat     = prctile(temp_img_grad_quadrat(:), 75); % max(temp_img_grad_quadrat(:)); % D/R
%             for k = 1:n_quadrat_x
%                 for l = 1:n_quadrat_y
%                     if temp_img_count_quadrat(k,l) == 0
%                         temp_img_dose_quadrat(k,l) = max_img_dose_quadrat;
%                         temp_img_peakdist_quadrat(k,l) = min_img_peakdist_quadrat;
%                         temp_img_grad_quadrat(k,l) = min_img_grad_quadrat;
%                     end
%                 end
%             end

            % Get relevant pixel and quadrat values
            temp_img_dose_quadrat     = temp_img_dose_quadrat .* temp_bw_quadrat;
            temp_img_peakdist_quadrat = temp_img_peakdist_quadrat .* temp_bw_quadrat;
            temp_img_grad_quadrat     = temp_img_grad_quadrat .* temp_bw_quadrat;

            PV_dose             = interp_img_dose(temp_bw);
            PV_dose_quadrat     = temp_img_dose_quadrat(temp_bw_quadrat); 
            PV_peakdist_quadrat = temp_img_peakdist_quadrat(temp_bw_quadrat);
            PV_grad_quadrat     = temp_img_grad_quadrat(temp_bw_quadrat);
            PV_count_quadrat    = temp_img_count_quadrat(temp_bw_quadrat);
            if any(strcmp(assay_pattern{i}, {'GRID_Stripes', 'GRID_Dots'}))
                PV_dose_peak_quadrat       = temp_img_dose_quadrat(temp_bw_peak_quadrat);
                PV_peakdist_peak_quadrat   = temp_img_peakdist_quadrat(temp_bw_peak_quadrat);
                PV_grad_peak_quadrat       = temp_img_grad_quadrat(temp_bw_peak_quadrat);
                PV_count_peak_quadrat      = temp_img_count_quadrat(temp_bw_peak_quadrat);
                PV_dose_valley_quadrat     = temp_img_dose_quadrat(temp_bw_valley_quadrat);
                PV_peakdist_valley_quadrat = temp_img_peakdist_quadrat(temp_bw_valley_quadrat);
                PV_grad_valley_quadrat     = temp_img_grad_quadrat(temp_bw_valley_quadrat);
                PV_count_valley_quadrat    = temp_img_count_quadrat(temp_bw_valley_quadrat);
            end

            % Store quadrat results
            assay.(assay_pattern{i}).img{j}.dose_avg             = mean(PV_dose);
            assay.(assay_pattern{i}).img{j}.N_colonies           = f(mean(PV_dose)) * N_colonies;
            assay.(assay_pattern{i}).img{j}.N_quadrat            = nnz(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.bw_dose              = bw_dose;
            assay.(assay_pattern{i}).img{j}.bw_dose_quadrat      = temp_bw_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.bw_quadrat           = temp_bw_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose             = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_peakdist         = interp_img_peakdist;
            assay.(assay_pattern{i}).img{j}.img_grad             = interp_img_grad;
            assay.(assay_pattern{i}).img{j}.img_dose_quadrat     = temp_img_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_peakdist_quadrat = temp_img_peakdist_quadrat;
            assay.(assay_pattern{i}).img{j}.img_grad_quadrat     = temp_img_grad_quadrat;
            assay.(assay_pattern{i}).img{j}.img_count_quadrat    = temp_img_count_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.dose_quadrat      = PV_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.peakdist_quadrat  = PV_peakdist_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.grad_quadrat      = PV_grad_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.count_quadrat     = PV_count_quadrat;
            if any(strcmp(assay_pattern{i}, {'GRID_Stripes', 'GRID_Dots'}))
                assay.(assay_pattern{i}).img{j}.N_quadrat_peak              = nnz(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.N_quadrat_valley            = nnz(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.bw_peak_quadrat             = temp_bw_peak_quadrat;
                assay.(assay_pattern{i}).img{j}.bw_valley_quadrat           = temp_bw_valley_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.dose_peak_quadrat        = PV_dose_peak_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.peakdist_peak_quadrat    = PV_peakdist_peak_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.grad_peak_quadrat        = PV_grad_peak_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.count_peak_quadrat       = PV_count_peak_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.dose_valley_quadrat      = PV_dose_valley_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.peakdist_valley_quadrat  = PV_peakdist_valley_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.grad_valley_quadrat      = PV_grad_valley_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.count_valley_quadrat     = PV_count_valley_quadrat;
            end
            
            % Plot quadrat images
            h = figure();
            set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
            tlo = tiledlayout(1,4);
            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_dose_quadrat, ...
                [0 dose_assay(j)+1], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_peakdist_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Distance to peak (cm)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_grad_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Dose gradient (cm/Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_count_quadrat, ...
                [0 max(temp_img_count_quadrat(:))], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Surviving colony count';
            hold on
            visboundaries(temp_bw_quadrat, 'Color', 'yellow', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)
            
%             % Plot image
%             h = figure();
%             set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%             imshow(assay.(assay_pattern{i}).img{j}.img_dose, ...
%                 [0 dose_assay(j)+1], 'colormap', jet(4096))
%             shading interp
%             c = colorbar;
%             c.Label.String = 'Dose (Gy)';
%             % title([assay.(assay_pattern{i}).img{j}.filename ': ' num2str(dose_assay(j)) ' Gy'])
%             hold on
%             for row = 1:dxdy_px(1):size(assay.(assay_pattern{i}).img{j}.img_dose,1)
%                 line([1, size(assay.(assay_pattern{i}).img{j}.img_dose,2)], ...
%                     [row, row], 'Color', [.7 .7 .7])
%             end
%             for col = 1:dxdy_px(2):size(assay.(assay_pattern{i}).img{j}.img_dose,2)
%                 line([col, col], [1, size(assay.(assay_pattern{i}).img{j}.img_dose,1)], ...
%                     'Color', [.7 .7 .7])
%             end
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask, ...
%                 'Color', 'yellow', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose, ...
%                 'Color', 'green', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_seg, ...
%                 'Color', 'red', 'LineWidth', 1)
%             plot(assay.(assay_pattern{i}).img{j}.y_centroids, ...
%                 assay.(assay_pattern{i}).img{j}.x_centroids, 'rx')
%             %             if ~isempty(xy_closestpeak)
%             %                 for ii = 1:length(stats)
%             %                     line([y_centroids(ii), xy_closestpeak(ii,1)], ...
%             %                         [x_centroids(ii), xy_closestpeak(ii,2)], ...
%             %                         'LineStyle', '-', 'LineWidth', 0.3, 'Color', 'r');
%             %                 end
%             %             end
%             hold off
%             set(gca, 'FontSize', 16)

%             figure();
%             imshow(assay.(assay_pattern{i}).img{j}.img_dose, ...
%                 [0 dose_assay(j)+1], 'colormap', jet(4096))
%             shading interp
%             c = colorbar;
%             c.Label.String = 'Dose (Gy)';
%             % title([assay.(assay_pattern{i}).img{j}.filename ': ' num2str(dose_assay(j)) ' Gy'])
%             hold on
%             for row = 1:dxdy_px(1):size(assay.(assay_pattern{i}).img{j}.img_dose,1)
%                 line([1, size(assay.(assay_pattern{i}).img{j}.img_dose,2)], ...
%                     [row, row], 'Color', [.7 .7 .7])
%             end
%             for col = 1:dxdy_px(2):size(assay.(assay_pattern{i}).img{j}.img_dose,2)
%                 line([col, col], [1, size(assay.(assay_pattern{i}).img{j}.img_dose,1)], ...
%                     'Color', [.7 .7 .7])
%             end
%             visboundaries(temp_bw_dose_peak, ...
%                 'Color', 'm', 'LineWidth', 1)
%             visboundaries(temp_bw_dose_valley, ...
%                 'Color', 'green', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_seg, ...
%                 'Color', 'red', 'LineWidth', 1)
%             plot(assay.(assay_pattern{i}).img{j}.y_centroids, ...
%                 assay.(assay_pattern{i}).img{j}.x_centroids, 'rx')
%             %             if ~isempty(xy_closestpeak)
%             %                 for ii = 1:length(stats)
%             %                     line([y_centroids(ii), xy_closestpeak(ii,1)], ...
%             %                         [x_centroids(ii), xy_closestpeak(ii,2)], ...
%             %                         'LineStyle', '-', 'LineWidth', 0.3, 'Color', 'r');
%             %                 end
%             %             end
%             hold off
%             set(gca, 'FontSize', 16)

        end
    end
end

%% Loop over controls first

N_colonies_ctrl           = [];
PE_ctrl                   = [];
avg_count_ctrl_quadrat    = [];

for i = 1:length(assay.Control.img)

    N_colonies_ctrl = [N_colonies_ctrl assay.Control.img{i}.N_colonies];
    PE_ctrl         = [PE_ctrl assay.Control.img{i}.PE];
    avg_count_ctrl_quadrat = [avg_count_ctrl_quadrat ...
        mean(assay.Control.img{i}.PV.count_quadrat)];

end

% [min(PE_ctrl) max(PE_ctrl)].*100
N_colonies_ctrl         = mean(N_colonies_ctrl);
PE_ctrl                 = mean(PE_ctrl);
avg_count_ctrl_quadrat  = mean(avg_count_ctrl_quadrat);

%% K-fold cross-validation for Poisson regression of quadrat data (training)

close all

path        = 'C:\Users\delmo\Desktop\New folder1\';

% Generate K-fold indices
K = 4;
cv_indices      = crossvalind('Kfold', 12, K);
cv_indices_ctrl = crossvalind('Kfold', 4, K);

% lgdstr  = {'Control', 'Dots', 'Dots - peak', 'Dots - valley', ...
%     'Stripes', 'Stripes - peak', 'Stripes - valley', 'Open', ...
%     'LQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
%     'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', ...
%     'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA)})', ...
%     'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR^2)})', ...
%     'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA^2)})'};

lgdstr  = {'Control', 'Dots', 'Dots - peak', 'Dots - valley', ...
    'Stripes', 'Stripes - peak', 'Stripes - valley', 'Open', ...
    'LQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
    'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', 'MLQ 95 % CI'};
color   = {'black', 'red', 'green', 'blue'};

significance    = 0.0001; % 0.000001 is too unstable for CI bands
assay_test      = {};
assay_train     = {};
alpha_LQ = []; alpha_MLQ_R = []; alpha_MLQ_A = []; alpha_MLQ_R2 = []; alpha_MLQ_A2 = [];
beta_LQ  = []; beta_MLQ_R  = []; beta_MLQ_A  = []; beta_MLQ_R2  = []; beta_MLQ_A2  = [];
delta_LQ   = []; delta_MLQ_R = []; delta_MLQ_A = []; delta_MLQ_R2 = []; delta_MLQ_A2 = [];

for k = 1:K

    test_indices        = (cv_indices == k);
    train_indices       = ~test_indices;
    test_indices        = find(test_indices);
    train_indices       = find(train_indices);
    test_indices_ctrl   = (cv_indices_ctrl == k);
    train_indices_ctrl  = ~test_indices_ctrl;
    test_indices_ctrl   = find(test_indices_ctrl);
    train_indices_ctrl  = find(train_indices_ctrl);

    for j = 1:length(assay_pattern)
        if ~strcmp(assay_pattern{j}, 'Control')
            for i = 1:length(test_indices)
                assay_test.(assay_pattern{j}).img{i} = ...
                    assay.(assay_pattern{j}).img{test_indices(i)};
            end
            for i = 1:length(train_indices)
                assay_train.(assay_pattern{j}).img{i} = ...
                    assay.(assay_pattern{j}).img{train_indices(i)};
            end
        else
            for i = 1:length(test_indices_ctrl)
                assay_test.(assay_pattern{j}).img{i} = ...
                    assay.(assay_pattern{j}).img{test_indices_ctrl(i)};
            end
            for i = 1:length(train_indices_ctrl)
                assay_train.(assay_pattern{j}).img{i} = ...
                    assay.(assay_pattern{j}).img{train_indices_ctrl(i)};
            end
        end
    end

    % Design variable matrix X
    X_D     = [];
    X_R     = [];
    X_G     = [];
    Y_CC    = [];
    Y_CTRL  = [];
    X_Open  = [];
    Y_Open  = [];

    % Structure quadrat dataset
    for j = 1:length(assay_pattern)

        dose_vec.(assay_pattern{j})      = [];
        peakdist_vec.(assay_pattern{j})  = [];
        grad_vec.(assay_pattern{j})      = [];
        count_vec.(assay_pattern{j})     = [];
        SF_vec.(assay_pattern{j})        = [];

        for i = 1:length(assay_train.(assay_pattern{j}).img)

            dose_quadrat_PV     = assay_train.(assay_pattern{j}).img{i}.PV.dose_quadrat;
            peakdist_quadrat_PV = assay_train.(assay_pattern{j}).img{i}.PV.peakdist_quadrat;
            grad_quadrat_PV     = assay_train.(assay_pattern{j}).img{i}.PV.grad_quadrat;
            count_quadrat_PV    = assay_train.(assay_pattern{j}).img{i}.PV.count_quadrat;

            % Compute surviving fraction (SF) per quadrat per irradiation
            % field configuration
            %             SF_quadrat = (count_quadrat .* f(dose_quadrat)) ./ ...
            %                 (avg_colonycount_ctrl_quadrat .* f(0));
            %             SF_quadrat = count_quadrat / avg_colonycount_ctrl_quadrat;

            dose_vec.(assay_pattern{j})     = [dose_vec.(assay_pattern{j}) dose_quadrat_PV'];
            peakdist_vec.(assay_pattern{j}) = [peakdist_vec.(assay_pattern{j}) peakdist_quadrat_PV'];
            grad_vec.(assay_pattern{j})     = [grad_vec.(assay_pattern{j}) grad_quadrat_PV'];
            count_vec.(assay_pattern{j})    = [count_vec.(assay_pattern{j}) count_quadrat_PV'];
            %             SF_vec.(assay_pattern{j})       = [SF_vec.(assay_pattern{j}) SF_quadrat];

        end

        % Design covariate matrices X
        if any(strcmp(assay_pattern{j}, {'Control', 'GRID_Dots', 'GRID_Stripes', 'Open'}))
            X_D     = [X_D      dose_vec.(assay_pattern{j})];
            X_R     = [X_R      peakdist_vec.(assay_pattern{j})];
            X_G     = [X_G      grad_vec.(assay_pattern{j})];
            Y_CC    = [Y_CC     count_vec.(assay_pattern{j})];
        end
        if strcmp(assay_pattern{j}, 'Control')
            Y_CTRL  = [Y_CTRL   count_vec.(assay_pattern{j})];
        end
        if any(strcmp(assay_pattern{j}, {'Control', 'Open'}))
            X_Open  = [X_Open   dose_vec.(assay_pattern{j})];
            Y_Open  = [Y_Open   count_vec.(assay_pattern{j})];
        end

    end

%     mdl_temp = fitglm(rand(length(Y_CTRL),1), Y_CTRL', 'y ~ x1', 'Distribution', 'poisson');
%     Y_CTRL   = random(mdl_temp, X_D(:));
    % Y_CTRL(Y_CTRL == 0) = 1;

%     [X_D, index_sort]       = sort(X_D);
%     X_R                     = X_R(index_sort);
%     X_A                     = X_A(index_sort);
%     Y_CC                    = Y_CC(index_sort);
%     [X_Open, index_sort]    = sort(X_Open);
%     Y_Open                  = Y_Open(index_sort);

    % Perform Poisson regression
    % alpha = 0.24 +- 0.02 Gy-1, beta = 0.019 +- 0.002 Gy-2
    tbl_LQ      = table(X_D(:), Y_CC(:), 'VariableNames', {'D', 'Y_CC'});
    tbl_MLQ_R   = table(X_D(:), X_R(:), Y_CC(:), 'VariableNames', {'D', 'R', 'Y_CC'});
    tbl_MLQ_G   = table(X_D(:), X_G(:), Y_CC(:), 'VariableNames', {'D', 'G', 'Y_CC'});
    tbl_MLQ_DG  = table(X_D(:), X_G(:), Y_CC(:), 'VariableNames', {'D', 'G', 'Y_CC'});
%     tbl_Open = table(X_Open', Y_Open', Y_CTRL', 'VariableNames', {'D_Open', 'Y_CC_Open', 'Y_CTRL'});

%     mdl_LQ = fitglm(tbl, 'Y_CC ~ 1 + D + D^2', 'Link', 'log', ...
%         'Distribution', 'poisson', 'Intercept', true);
%     mdl_MLQ_R = fitglm(tbl, 'Y_CC ~ 1 + D + D^2 + R', 'Link', 'log', ...
%         'Distribution', 'poisson', 'Intercept', true);
%     mdl_MLQ_G = fitglm(tbl, 'Y_CC ~ 1 + D + D^2 + G', 'Link', 'log', ...
%         'Distribution', 'poisson', 'Intercept', true);
%     mdl_MLQ_DG = fitglm(tbl, 'Y_CC ~ 1 + D + D^2 + D:G', 'Link', 'log', ...
%         'Distribution', 'poisson', 'Intercept', true);
%     mdl_Open = fitglm(tbl_Open, 'Y_CC_Open ~ 1 + D_Open + D_Open^2', ...
%         'Link', 'log', 'Distribution', 'poisson', 'Intercept', true);

    mdl_LQ = fitglm(tbl_LQ, 'Y_CC ~ D + D^2', ...
        'Link', 'log', 'Distribution', 'poisson'); % 'Intercept', false
    mdl_MLQ_R = fitglm(tbl_MLQ_R, 'Y_CC ~ D + D^2 + R', ...
        'Link', 'log', 'Distribution', 'poisson'); % 'Intercept', false
    mdl_MLQ_G = fitglm(tbl_MLQ_G, 'Y_CC ~ D + D^2 + G', ...
        'Link', 'log', 'Distribution', 'poisson');
    mdl_MLQ_DG = fitglm(tbl_MLQ_DG, 'Y_CC ~ D + D^2 + D:G', ...
        'Link', 'log', 'Distribution', 'poisson');
%     mdl_Open = fitglm(tbl_Open, 'Y_CC_Open ~ D_Open + D_Open^2', ...
%         'Distribution', 'poisson', 'Offset', 'Y_CTRL');

    [mdl_LQ.ModelCriterion.AIC mdl_LQ.ModelCriterion.BIC; ...
        mdl_MLQ_R.ModelCriterion.AIC mdl_MLQ_R.ModelCriterion.BIC; ...
        mdl_MLQ_G.ModelCriterion.AIC mdl_MLQ_G.ModelCriterion.BIC; ...
        mdl_MLQ_DG.ModelCriterion.AIC mdl_MLQ_DG.ModelCriterion.BIC]

%     G = mdl_LQ.Deviance - mdl_MLQ_R.Deviance;
%     p = 1 - chi2cdf(G, 1)

    % Create data points for prediction
    D_fit            = linspace(min(X_D(:)), max(X_D(:)), 101);
    R_fit            = linspace(min(X_R(:)), max(X_R(:)), 101);
    G_fit            = linspace(min(X_G(:)), max(X_G(:)), 101);
    % CTRL_fit         = linspace(min(Y_CTRL(:)), max(Y_CTRL(:)), 101);
    [D_mesh1, R_mesh] = meshgrid(D_fit, R_fit);  
    [D_mesh2, G_mesh] = meshgrid(D_fit, G_fit);

    isequal(D_mesh1, D_mesh2)

    ypred_LQ     = predict(mdl_LQ, [D_mesh(:)]);
    ypred_MLQ_R  = predict(mdl_MLQ_R, [D_mesh(:), R_mesh(:)]);
    ypred_MLQ_G  = predict(mdl_MLQ_G, [D_mesh(:), G_mesh(:)]);
    ypred_MLQ_DG = predict(mdl_MLQ_DG, [D_mesh(:), G_mesh(:)]);

    % Plot prediction surface
    plot3Dsurf(X_D, X_R, Y_CC, D_mesh, R_mesh, ...
        reshape((ypred_LQ),101,101), reshape((ypred_MLQ_R),101,101), ...
        'Distance to peak (cm)', ...
        'MLQ fit $(\mu_0 \exp(b_0 -\alpha D - \beta D^2 + \delta R))$', ...
        fullfile(path, sprintf('SurfPlot_LQvsMLQ_R_%s.png', num2str(k))))
    plot3Dsurf(X_D, X_G, Y_CC, D_mesh, G_mesh, ypred_LQ, ypred_MLQ_G, ...
        'Dose gradient (cm/Gy)', ...
        'MLQ fit $(\mu_0 \exp(-\alpha D - \beta D^2 + \delta \bar{\dot{D}}))$', ...
        fullfile(path, sprintf('SurfPlot_LQvsMLQ_G_%s.png', num2str(k))))
    
%     plot3Dsurf(X_D, X_R, Y_CC, D_mesh, G_mesh, ypred_LQ, ypred_MLQ_DG, ...
%         'Dose gradient (cm/Gy)', ...
%         'MLQ fit $(\mu_0 \exp(-\alpha D -\beta D^2 + \delta D \bar{\dot{D}}))$', ...
%         fullfile(path, sprintf('SurfPlot_LQvsMLQ_DG_%s.png', num2str(k))))

end
