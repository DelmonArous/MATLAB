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

    % Compute gradient image (cm/Gy)
%     img_grad                     = img_peakdist ./ temp_img_dose;
%     img_grad(temp_bw_dose == 0)  = 0;
%     img_grad(img_grad <= 0)      = 0;
%     img_grad(isnan(img_grad))    = 0;
%     img_grad(isinf(img_grad))    = 0;

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
            bw_assay_seg = imtranslate(bw_assay_seg, ...
                translation_vec.(assay_pattern{i}){j});
            bw_flask = imtranslate(bw_flask, ...
                translation_vec.(assay_pattern{i}){j});
            %                         figure(); imshow(bw_assay_seg, [])
            %                         figure(); imshow(bw_flask, [])
            %         elseif contains(pattern_name, "Open") || ...
            %                 (contains(pattern_name, ["Dots", "Stripes"]) && j <= 4)
            %             filename
            %             figure(); imshow(bw_flask, [])
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

% %         % 1D survival analysis
% %         x_centroid_1D   = sort(x_centroids);
% % 
% %         colony_per_pxrow = [];
% %         for k = 1:size(bw_flask, 1)
% %             colony_per_pxrow = [colony_per_pxrow; ...
% %                 length(x_centroid_1D(x_centroid_1D == k))];
% %         end
% %         row_weights      = sum(bw_flask, 2);
% %         row_weights      = row_weights ./ max(row_weights);
% %         colony_per_pxrow = colony_per_pxrow ./ row_weights;
% %         colony_per_pxrow(isinf(colony_per_pxrow)) = 0;
% %         colony_per_pxrow(isnan(colony_per_pxrow)) = 0;
% % 
% %         n_band              = ceil(size(bw_flask,1)./dx_px);
% %         dummy_zeros         = zeros(n_band*dx_px - size(bw_flask, 1), 1);
% %         colony_per_pxrow    = cat(1, colony_per_pxrow, dummy_zeros);
% %         colony_count_band   = sum(reshape(colony_per_pxrow, ...
% %             length(colony_per_pxrow)/n_band, n_band));
% %         position_band       = linspace(0, size(bw_flask, 1)*px_size_assay, n_band);
% % 
% %         % Store band results
% %         assay.(assay_pattern{i}).img{j}.position_band     = position_band;
% %         assay.(assay_pattern{i}).img{j}.colony_count_band = colony_count_band;

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

        % Compute peak distance for each centroid
%         if any(strcmp(assay_pattern{i}, {'GRID_Dots', 'GRID_Stripes'}))
%             [img_peakdist, xy_closestpeak] = EstimatePeakDistanceMap( ...
%                 bw_centroids, centroids, peak_xycoord.(assay_pattern{i}));
%         else
%             img_peakdist   = zeros(size(bw_assay_seg,1), size(bw_assay_seg,2));
%             xy_closestpeak = [];
%         end
%         img_peakdist = px_size_assay .* img_peakdist;
        
        % Define peak area image
%         img_grad = zeros(size(bw_assay_seg,1), size(bw_assay_seg,2));

        % Dividing the segmentation and flask mask into
        % (n_quadrat_y x n_quadrat_x) quadrats
        bw_flask_quadrat       = mat2quadrat(bw_flask, dxdy_px);
        bw_assay_seg_quadrat   = mat2quadrat(bw_assay_seg, dxdy_px);
        bw_centroids_quadrat   = mat2quadrat(bw_centroids, dxdy_px);

        temp_bw_flask_quadrat      = false(n_quadrat_x, n_quadrat_y);
        temp_img_count_quadrat     = zeros(n_quadrat_x, n_quadrat_y);

%         img_peakdist_quadrat   = mat2quadrat(img_peakdist, dxdy_px);
%         img_grad_quadrat       = mat2quadrat(img_grad, dxdy_px);
%         temp_img_peakdist_quadrat  = zeros(n_quadrat_x, n_quadrat_y);
%         temp_img_grad_quadrat      = zeros(n_quadrat_x, n_quadrat_y);

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

%                 if ~isempty(x_temp)
%                     centroid_quadrat{k,l}.peak_distance = ...
%                         sum(img_peakdist_quadrat{k,l}(:)) / length(x_temp);
% %                     centroid_quadrat{k,l}.grad          = ...
% %                         sum(img_grad_quadrat{k,l}(:)) / length(x_temp);
%                 else
%                     centroid_quadrat{k,l}.peak_distance = 0.0;
% %                     centroid_quadrat{k,l}.grad          = 0.0;
%                 end
%                 temp_img_peakdist_quadrat(k,l) = centroid_quadrat{k,l}.peak_distance;
%                 temp_img_grad_quadrat(k,l)      = centroid_quadrat{k,l}.grad;

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
%         assay.(assay_pattern{i}).img{j}.img_peakdist_quadrat = temp_img_peakdist_quadrat;
%         assay.(assay_pattern{i}).img{j}.img_grad_quadrat     = temp_img_grad_quadrat;

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
%             interp_img_grad = interp2(X, Y, temp_img_grad, ...
%                 X2, Y2, 'spline');
%             interp_img_grad = cat(1, interp_img_grad, dummy_zeros_x);
%             interp_img_grad = cat(2, interp_img_grad, dummy_zeros_y);

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
%                 interp_img_grad = imwarp(interp_img_grad, tform, ...
%                     'OutputView', Rfixed);

            elseif j > 1 && strcmp(assay_pattern{i}, 'Open')

                bw_dose = imwarp(bw_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_dose = imwarp(interp_img_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_peakdist = imwarp(interp_img_peakdist, tform, ...
                    'OutputView', Rfixed);
%                 interp_img_grad = imwarp(interp_img_grad, tform, ...
%                     'OutputView', Rfixed);

            elseif any(strcmp(assay_pattern{i}, {'GRID_Stripes', 'GRID_Dots'}))

                bw_dose = imtranslate(bw_dose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_dose = imtranslate(interp_img_dose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_peakdist = imtranslate(interp_img_peakdist, ...
                    translation_dose_vec.(assay_pattern{i}));
%                 interp_img_grad = imtranslate(interp_img_grad, ...
%                     translation_dose_vec.(assay_pattern{i}));

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
            interp_img_grad = interp_img_peakdist ./ interp_img_dose;
            temp_img_grad_quadrat = temp_img_peakdist_quadrat ./ temp_img_dose_quadrat;
            temp_img_grad_quadrat(temp_bw_dose_quadrat == 0)    = 0;
            temp_img_grad_quadrat(temp_img_grad_quadrat <= 0)   = 0;
            temp_img_grad_quadrat(isnan(temp_img_grad_quadrat)) = 0;
            temp_img_grad_quadrat(isinf(temp_img_grad_quadrat)) = 0;

            temp_bw                = bw_flask & bw_dose;
            temp_bw_peak           = temp_bw_dose_peak & temp_bw;
            temp_bw_valley         = temp_bw_dose_valley & temp_bw;
            temp_bw_quadrat        = temp_bw_flask_quadrat & temp_bw_dose_quadrat;
            temp_bw_peak_quadrat   = temp_bw_dose_peak_quadrat & temp_bw_quadrat;
            temp_bw_valley_quadrat = temp_bw_dose_valley_quadrat & temp_bw_quadrat;
%             temp_bw_quadrat = temp_bw_dose_peak_quadrat | temp_bw_dose_valley_quadrat; % NEW!!!

            % Get relevant pixel and quadrat values
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
            assay.(assay_pattern{i}).img{j}.dose_avg            = mean(PV_dose);
            assay.(assay_pattern{i}).img{j}.N_colonies          = f(mean(PV_dose)) * N_colonies;
            assay.(assay_pattern{i}).img{j}.N_quadrat           = nnz(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.bw_dose             = bw_dose;
            assay.(assay_pattern{i}).img{j}.bw_dose_quadrat     = temp_bw_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.bw_quadrat          = temp_bw_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose            = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_peakdist        = interp_img_peakdist;
            assay.(assay_pattern{i}).img{j}.img_grad            = interp_img_grad;
            assay.(assay_pattern{i}).img{j}.img_dose_quadrat    = temp_img_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_grad_quadrat    = temp_img_grad_quadrat;
            assay.(assay_pattern{i}).img{j}.img_count_quadrat   = temp_img_count_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.dose_quadrat     = PV_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.peakdist_quadrat = PV_peakdist_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.grad_quadrat     = PV_grad_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.count_quadrat    = PV_count_quadrat;
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

            h = figure();
            set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
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

            h = figure();
            set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
            imshow(assay.(assay_pattern{i}).img{j}.img_count_quadrat, ...
                [0 max(temp_img_count_quadrat(:))], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Surviving colony count';
            hold on
            visboundaries(temp_bw_quadrat, 'Color', 'yellow', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            figure();
            imshow(assay.(assay_pattern{i}).img{j}.img_dose, ...
                [0 dose_assay(j)+1], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Dose (Gy)';
            % title([assay.(assay_pattern{i}).img{j}.filename ': ' num2str(dose_assay(j)) ' Gy'])
            hold on
            for row = 1:dxdy_px(1):size(assay.(assay_pattern{i}).img{j}.img_dose,1)
                line([1, size(assay.(assay_pattern{i}).img{j}.img_dose,2)], ...
                    [row, row], 'Color', [.7 .7 .7])
            end
            for col = 1:dxdy_px(2):size(assay.(assay_pattern{i}).img{j}.img_dose,2)
                line([col, col], [1, size(assay.(assay_pattern{i}).img{j}.img_dose,1)], ...
                    'Color', [.7 .7 .7])
            end
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_seg, ...
                'Color', 'red', 'LineWidth', 1)
            plot(assay.(assay_pattern{i}).img{j}.y_centroids, ...
                assay.(assay_pattern{i}).img{j}.x_centroids, 'rx')
            %             if ~isempty(xy_closestpeak)
            %                 for ii = 1:length(stats)
            %                     line([y_centroids(ii), xy_closestpeak(ii,1)], ...
            %                         [x_centroids(ii), xy_closestpeak(ii,2)], ...
            %                         'LineStyle', '-', 'LineWidth', 0.3, 'Color', 'r');
            %                 end
            %             end
            hold off
            set(gca, 'FontSize', 16)

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

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_quadrat_vec     = [];
% dose_quadrat_vec  = [];
% count_quadrat_vec = []; count_quadrat_sd_vec = [];
% 
% N_ctrl_quadrat_vec     = [];
% dose_ctrl_quadrat_vec  = [];
% count_ctrl_quadrat_vec = []; count_ctrl_quadrat_sd_vec = [];
% 
% N_peak_dots_quadrat_vec     = [];
% dose_peak_dots_quadrat_vec  = [];
% count_peak_dots_quadrat_vec = []; count_peak_dots_quadrat_sd_vec = [];
% 
% N_valley_dots_quadrat_vec     = [];
% dose_valley_dots_quadrat_vec  = [];
% count_valley_dots_quadrat_vec = []; count_valley_dots_quadrat_sd_vec = [];
% 
% N_peak_stripes_quadrat_vec     = [];
% dose_peak_stripes_quadrat_vec  = [];
% count_peak_stripes_quadrat_vec = []; count_peak_stripes_quadrat_sd_vec = [];
% 
% N_valley_stripes_quadrat_vec     = [];
% dose_valley_stripes_quadrat_vec  = [];
% count_valley_stripes_quadrat_vec = []; count_valley_stripes_quadrat_sd_vec = [];
% 
% for j = 1:length(assay_pattern)
% 
%     for i = 1:length(assay.(assay_pattern{j}).img)
% 
%         if strcmp(assay_pattern{j}, {'Control'})
%             N_ctrl_quadrat_vec = [N_ctrl_quadrat_vec ...
%                 assay.(assay_pattern{j}).img{i}.N_quadrat];
%             dose_ctrl_quadrat_vec = [dose_ctrl_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_quadrat)];
%             count_ctrl_quadrat_vec = [count_ctrl_quadrat_vec ...
%                 sum(assay.(assay_pattern{j}).img{i}.PV.count_quadrat)];
%         elseif strcmp(assay_pattern{j}, {'Open'})
%             N_quadrat_vec = [N_quadrat_vec ...
%                 assay.(assay_pattern{j}).img{i}.N_quadrat];
%             dose_quadrat_vec = [dose_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_quadrat)];
%             count_quadrat_vec = [count_quadrat_vec ...
%                 sum(assay.(assay_pattern{j}).img{i}.PV.count_quadrat)];
%         elseif strcmp(assay_pattern{j}, {'GRID_Stripes'})
%             N_peak_stripes_quadrat_vec = [N_peak_stripes_quadrat_vec ...
%                 assay.(assay_pattern{j}).img{i}.N_quadrat_peak];
%             dose_peak_stripes_quadrat_vec = [dose_peak_stripes_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_peak_quadrat)];
%             count_peak_stripes_quadrat_vec = [count_peak_stripes_quadrat_vec ...
%                 sum(assay.(assay_pattern{j}).img{i}.PV.count_peak_quadrat)];
%             N_valley_stripes_quadrat_vec = [N_valley_stripes_quadrat_vec ...
%                 assay.(assay_pattern{j}).img{i}.N_quadrat_valley];
%             dose_valley_stripes_quadrat_vec = [dose_valley_stripes_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_valley_quadrat)];
%             count_valley_stripes_quadrat_vec = [count_valley_stripes_quadrat_vec ...
%                 sum(assay.(assay_pattern{j}).img{i}.PV.count_valley_quadrat)];
%         elseif strcmp(assay_pattern{j}, {'GRID_Dots'})
%             N_peak_dots_quadrat_vec = [N_peak_dots_quadrat_vec ...
%                 assay.(assay_pattern{j}).img{i}.N_quadrat_peak];
%             dose_peak_dots_quadrat_vec = [dose_peak_dots_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_peak_quadrat)];
%             count_peak_dots_quadrat_vec = [count_peak_dots_quadrat_vec ...
%                 sum(assay.(assay_pattern{j}).img{i}.PV.count_peak_quadrat)];
%             N_valley_dots_quadrat_vec = [N_valley_dots_quadrat_vec ...
%                 assay.(assay_pattern{j}).img{i}.N_quadrat_valley];
%             dose_valley_dots_quadrat_vec = [dose_valley_dots_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_valley_quadrat)];
%             count_valley_dots_quadrat_vec = [count_valley_dots_quadrat_vec ...
%                 sum(assay.(assay_pattern{j}).img{i}.PV.count_valley_quadrat)];
%         end
% 
%     end
% 
%     if strcmp(assay_pattern{j}, {'Control'})
%         count_ctrl_quadrat_vec = count_ctrl_quadrat_vec./N_ctrl_quadrat_vec/avg_count_ctrl_quadrat;
%     elseif strcmp(assay_pattern{j}, {'Open'})
%         count_quadrat_vec = count_quadrat_vec./N_quadrat_vec/avg_count_ctrl_quadrat;
%     elseif strcmp(assay_pattern{j}, {'GRID_Stripes'})
%         count_peak_stripes_quadrat_vec = count_peak_stripes_quadrat_vec./N_peak_stripes_quadrat_vec/avg_count_ctrl_quadrat;
%         count_valley_stripes_quadrat_vec = count_valley_stripes_quadrat_vec./N_valley_stripes_quadrat_vec/avg_count_ctrl_quadrat;
%     elseif strcmp(assay_pattern{j}, {'GRID_Dots'})
%         count_peak_dots_quadrat_vec = count_peak_dots_quadrat_vec./N_peak_dots_quadrat_vec/avg_count_ctrl_quadrat;
%         count_valley_dots_quadrat_vec = count_valley_dots_quadrat_vec./N_valley_dots_quadrat_vec/avg_count_ctrl_quadrat;
%     end
% 
% end
% 
% % avg_count_ctrl_quadrat
% D = [dose_ctrl_quadrat_vec dose_quadrat_vec ...
%      dose_peak_stripes_quadrat_vec dose_valley_stripes_quadrat_vec ...
%      dose_peak_dots_quadrat_vec dose_valley_dots_quadrat_vec]';
% CC = [count_ctrl_quadrat_vec count_quadrat_vec ...
%     count_peak_stripes_quadrat_vec count_valley_stripes_quadrat_vec ...
%     count_peak_dots_quadrat_vec count_valley_dots_quadrat_vec]';
% % SF(isinf(SF)) = 0;
% 
% X = [D, D.^2];
% modelfun = @(b,x) b(1)*exp(b(2)*x(:,1) + b(3)*x(:,2));
% beta0 = [1 -0.1 -0.1];
% mdl = fitnlm(X, CC, modelfun, beta0);
% 
% % Regression
% % tbl = table(D, CC, 'VariableNames', {'D', 'CC'});
% % mdl = fitlm(tbl, 'CC ~ D + D^2', 'Intercept', true);
% a       = mdl.Coefficients.Estimate(1);
% alpha   = mdl.Coefficients.Estimate(2);
% beta    = mdl.Coefficients.Estimate(3);
% D_fit   = linspace(0, 10, 1001);
% CC_fit  = a.*exp(alpha.*D_fit + beta.*(D_fit.^2));
% 
% % Plot
% 
% h_vec   = [];
% lgdstr  = {'Control', 'Open', 'Stripes, peaks', 'Stripes, valleys', ...
%     'Dots, peaks', 'Dots, valleys'};
% color   = {'black', 'green', 'blue', 'red'};
% 
% figure(242);
% hold on
% h_ctrl = plot(dose_ctrl_quadrat_vec, count_ctrl_quadrat_vec, ...
%     'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{1});
% h_open = plot(dose_quadrat_vec, count_quadrat_vec, 'o', 'MarkerSize', 6, ...
%     'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{4});
% h_p_dots = plot(dose_peak_dots_quadrat_vec, count_peak_dots_quadrat_vec, ...
%     '^', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{2});
% h_v_dots = plot(dose_valley_dots_quadrat_vec, count_valley_dots_quadrat_vec, ...
%     'v', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{2});
% h_p_stripes = plot(dose_peak_stripes_quadrat_vec, count_peak_stripes_quadrat_vec, ...
%     '^', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{3});
% h_v_stripes = plot(dose_valley_stripes_quadrat_vec, count_valley_stripes_quadrat_vec, ...
%     'v', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{3});
% h_fit = plot(D_fit, CC_fit, 'r-');
% hold off
% h_vec = [h_ctrl h_open h_p_stripes h_v_stripes h_p_dots h_v_dots h_fit];
% legend(h_vec, [lgdstr, 'LQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2)})'], ...
%     'Location', 'NorthEast')
% xlabel('Average quadrat dose (Gy)')
% ylabel('Surviving fraction (per quadrat)')
% grid on
% set(gca, 'FontSize', 18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% temp_dose_vec     = [];
% temp_count_vec    = [];
% temp_count_sd_vec = [];
% 
% temp_peak_dose_vec     = [];
% temp_peak_count_vec    = [];
% temp_peak_count_sd_vec = [];
% 
% temp_valley_dose_vec     = [];
% temp_valley_count_vec    = [];
% temp_valley_count_sd_vec = [];
% 
% for j = 1:length(assay_pattern)
% 
%     N_quadrat_vec     = [];
%     dose_quadrat_vec  = [];
%     count_quadrat_vec = []; count_quadrat_sd_vec = [];
%     SF_quadrat_vec    = [];
% 
%     N_peak_quadrat_vec     = [];
%     dose_peak_quadrat_vec  = [];
%     count_peak_quadrat_vec = []; count_peak_quadrat_sd_vec = [];
%     SF_peak_quadrat_vec    = [];
% 
%     N_valley_quadrat_vec     = [];
%     dose_valley_quadrat_vec  = [];
%     count_valley_quadrat_vec = []; count_valley_quadrat_sd_vec = [];
%     SF_valley_quadrat_vec    = [];
% 
%     for i = 1:length(assay.(assay_pattern{j}).img)
%         N_quadrat_vec = [N_quadrat_vec ...
%             length(assay.(assay_pattern{j}).img{i}.PV.dose_quadrat)]; 
%         dose_quadrat_vec = [dose_quadrat_vec ...
%             mean(assay.(assay_pattern{j}).img{i}.PV.dose_quadrat)];
%         count_quadrat_vec = [count_quadrat_vec ...
%             mean(assay.(assay_pattern{j}).img{i}.PV.count_quadrat)];
%         count_quadrat_sd_vec = [count_quadrat_sd_vec ...
%             std(assay.(assay_pattern{j}).img{i}.PV.count_quadrat)];
% 
%         if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
%             N_peak_quadrat_vec = [N_peak_quadrat_vec ...
%                 length(assay.(assay_pattern{j}).img{i}.PV.dose_peak_quadrat)];
%             dose_peak_quadrat_vec = [dose_peak_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_peak_quadrat)];
%             count_peak_quadrat_vec = [count_peak_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.count_peak_quadrat)];
%             count_peak_quadrat_sd_vec = [count_peak_quadrat_sd_vec ...
%                 std(assay.(assay_pattern{j}).img{i}.PV.count_peak_quadrat)];
%             N_valley_quadrat_vec = [N_valley_quadrat_vec ...
%                 length(assay.(assay_pattern{j}).img{i}.PV.dose_valley_quadrat)];
%             dose_valley_quadrat_vec = [dose_valley_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.dose_valley_quadrat)];
%             count_valley_quadrat_vec = [count_valley_quadrat_vec ...
%                 mean(assay.(assay_pattern{j}).img{i}.PV.count_valley_quadrat)];
%             count_valley_quadrat_sd_vec = [count_valley_quadrat_sd_vec ...
%                 std(assay.(assay_pattern{j}).img{i}.PV.count_valley_quadrat)];
%         end
% 
%     end
% 
%     if strcmp(assay_pattern{j}, 'Control')
%         temp_dose_vec = [temp_dose_vec ...
%             sum(N_quadrat_vec.*dose_quadrat_vec)/sum(N_quadrat_vec)];
%         temp_count_vec = [temp_count_vec ...
%             sum(N_quadrat_vec.*count_quadrat_vec)/sum(N_quadrat_vec)];
%         %         temp_count_sd_vec = [temp_count_sd_vec ...
%         %             sqrt(sum((N_quadrat_vec-1).*count_quadrat_sd_vec.^2) / ...
%         %             (sum(N_quadrat_vec)-1))];
%         temp_count_sd_vec = [temp_count_sd_vec ...
%             std(count_quadrat_vec)];
%     elseif strcmp(assay_pattern{j}, 'Open')
%         for k = 1:3
%             temp_dose_vec = [temp_dose_vec ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)) .* ...
%                 dose_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)))];
%             temp_count_vec = [temp_count_vec ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)).* ...
%                 count_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)))];
%             %             temp_count_sd_vec = [temp_count_sd_vec ...
%             %                 sqrt(sum((N_quadrat_vec((4*k-3):(4*k))-1).* ...
%             %                 count_quadrat_sd_vec((4*k-3):(4*k)).^2) / ...
%             %                 (sum(N_quadrat_vec((4*k-3):(4*k)))-1))];
%             temp_count_sd_vec = [temp_count_sd_vec ...
%                 std(count_quadrat_vec((4*k-3):(4*k)))];
%         end
%     else
%         for k = 1:3
%             temp_dose_vec = [temp_dose_vec ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)) .* ...
%                 dose_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)))];
%             temp_count_vec = [temp_count_vec ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)).* ...
%                 count_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_quadrat_vec((4*k-3):(4*k)))];
%             temp_count_sd_vec = [temp_count_sd_vec ...
%                 std(count_quadrat_vec((4*k-3):(4*k)))];
%             temp_peak_dose_vec = [temp_peak_dose_vec ...
%                 sum(N_peak_quadrat_vec((4*k-3):(4*k)) .* ...
%                 dose_peak_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_peak_quadrat_vec((4*k-3):(4*k)))];
%             temp_peak_count_vec = [temp_peak_count_vec ...
%                 sum(N_peak_quadrat_vec((4*k-3):(4*k)).* ...
%                 count_peak_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_peak_quadrat_vec((4*k-3):(4*k)))];
%             temp_peak_count_sd_vec = [temp_peak_count_sd_vec ...
%                 std(count_peak_quadrat_vec((4*k-3):(4*k)))];
%             temp_valley_dose_vec = [temp_valley_dose_vec ...
%                 sum(N_valley_quadrat_vec((4*k-3):(4*k)) .* ...
%                 dose_valley_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_valley_quadrat_vec((4*k-3):(4*k)))];
%             temp_valley_count_vec = [temp_valley_count_vec ...
%                 sum(N_valley_quadrat_vec((4*k-3):(4*k)).* ...
%                 count_valley_quadrat_vec((4*k-3):(4*k))) / ...
%                 sum(N_valley_quadrat_vec((4*k-3):(4*k)))];
%             temp_valley_count_sd_vec = [temp_valley_count_sd_vec ...
%                 std(count_valley_quadrat_vec((4*k-3):(4*k)))];
%         end
% 
%     end
% 
% end

% h_vec   = [];
% lgdstr  = {'Control', 'Dots', 'Stripes', 'Open'};
% color   = {'black', 'blue', 'green', 'red'};
% facecolor = {[.7 .7 .7], [0.3010 0.7450 0.9330], [0.78 1.0 0.46], [1.0 0.8 0.8]};
% temp_dose_vec(5) = 0.8;
% x = categorical(string(round(temp_dose_vec,1)));
% x = reordercats(x, string(x));
% 
% figure(400);
% hold on 
% h = bar(x(1), temp_count_vec(1), color{1}, 'LineWidth', 2, ...
%     'FaceColor', facecolor{1}, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x(1), temp_count_vec(1), temp_count_sd_vec(1), ...
%     'Color', [0 0 0], 'LineStyle', 'none');
% h_vec = [h_vec h];
% for k = 1:3
%     h = bar(x((3*k-1):(3*k+1)), temp_count_vec((3*k-1):(3*k+1)), ...
%         color{k+1}, 'LineWidth', 2, 'FaceColor', facecolor{k+1}, ...
%         'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
%     errorbar(x((3*k-1):(3*k+1)), temp_count_vec((3*k-1):(3*k+1)), ...
%         temp_count_sd_vec((3*k-1):(3*k+1)), 'Color', [0 0 0], ...
%         'LineStyle', 'none')
%     h_vec = [h_vec h];
% %     lgdstr{k} = [lgdstr{k} '($N=$' num2str() ')'];
% end
% hold off
% legend(h_vec, lgdstr, 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% xlabel('Average quadrat dose (Gy)')
% ylabel('Surviving colony count (per quadrat)')
% grid on
% set(gca, 'FontSize', 18)

% h_vec   = [];
% lgdstr  = {'Peak, Dots', 'Valley, Dots', 'Peak, Stripes', 'Valley, Stripes'};
% color   = {'blue', 'green'};
% facecolor = {[0.78 1.0 0.46], [1.0 0.8 0.8]};
% hatchtype = {'speckle', 'single'};
% % temp_dose_vec = [temp_peak_dose_vec temp_valley_dose_vec];
% % temp_count_vec = [temp_peak_count_vec temp_valley_count_vec];
% % temp_count_sd_vec = [temp_peak_count_sd_vec temp_valley_count_sd_vec];
% x_peak = categorical(string(round(temp_peak_dose_vec,1)));
% x_peak = reordercats(x_peak, string(x_peak));
% x_valley = categorical(string(round(temp_valley_dose_vec,1)));
% x_valley = reordercats(x_valley, string(x_valley));
% f = figure(401);
% hold on
% for k = 1:2
%     h_peak = bar(x_peak((3*k-2):(3*k)), temp_peak_count_vec((3*k-2):(3*k)), ...
%         color{k}, 'LineWidth', 2, 'FaceColor', facecolor{1}, ...
%         'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
%     errorbar(x_peak((3*k-2):(3*k)), temp_peak_count_vec((3*k-2):(3*k)), ...
%         temp_peak_count_sd_vec((3*k-2):(3*k)), 'Color', [0 0 0], ...
%         'LineStyle', 'none')
%     h_valley = bar(x_valley((3*k-2):(3*k)), temp_valley_count_vec((3*k-2):(3*k)), ...
%         color{k}, 'LineWidth', 2, 'FaceColor', facecolor{2}, ...
%         'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
%     errorbar(x_valley((3*k-2):(3*k)), temp_valley_count_vec((3*k-2):(3*k)), ...
%         temp_valley_count_sd_vec((3*k-2):(3*k)), 'Color', [0 0 0], ...
%         'LineStyle', 'none')
%     h_vec = [h_vec h_peak h_valley];
% %     lgdstr{k} = [lgdstr{k} '($N=$' num2str() ')']; % Legg til antall
% %     kvadranter her for 
% end
% hold off
% legend(h_vec, lgdstr, 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% xlabel('Average quadrat dose (Gy)')
% ylabel('Surviving colony count (per quadrat)')
% ylim([0 Inf])
% grid on
% set(gca, 'FontSize', 18)

% lgdstr  = {'Control', 'Open', 'Valley, Dots', 'Valley, Stripes'};
% ttemp_dose_vec = [temp_dose_vec(1) temp_valley_dose_vec temp_dose_vec(8)];
% ttemp_count_vec = [temp_count_vec(1) temp_valley_count_vec temp_count_vec(8)];
% ttemp_count_sd_vec = [temp_count_sd_vec(1) temp_valley_count_sd_vec temp_count_sd_vec(8)];
% x = categorical(string(round(ttemp_dose_vec, 1)));
% 
% figure(401);
% hold on
% h_ctrl = bar(x(1), ttemp_count_vec(1), 'k', 'LineWidth', 2, ...
%     'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x(1), ttemp_count_vec(1), ttemp_count_sd_vec(1), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% h_open = bar(x(end), ttemp_count_vec(end), 'r', 'LineWidth', 2, ...
%     'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x(end), ttemp_count_vec(end), ttemp_count_sd_vec(end), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% h_dots_valley = bar(x(2:4), ttemp_count_vec(2:4), 'g', ...
%     'LineWidth', 2, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x(2:4), ttemp_count_vec(2:4), ttemp_count_sd_vec(2:4), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% h_stripes_valley = bar(x(5:7), ttemp_count_vec(5:7), 'b', ...
%     'LineWidth', 2, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x(5:7), ttemp_count_vec(5:7), ttemp_count_sd_vec(5:7), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% %     lgdstr{k} = [lgdstr{k} '($N=$' num2str() ')']; % Legg til antall kvadranter her for
% hold off
% legend([h_ctrl h_open h_dots_valley h_stripes_valley], lgdstr, ...
%     'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% xlabel('Average quadrat dose (Gy)')
% ylabel('Surviving colony count (per quadrat)')
% ylim([0 Inf])
% grid on
% set(gca, 'FontSize', 18)
% 
% 
% lgdstr  = {'Open', 'Peak, Dots', 'Peak, Stripes'};
% ttemp_dose_vec = [temp_peak_dose_vec temp_dose_vec(8:10)];
% ttemp_count_vec = [temp_peak_count_vec temp_count_vec(8:10)];
% ttemp_count_sd_vec = [temp_peak_count_sd_vec temp_count_sd_vec(8:10)];
% [ttemp_dose_vec, ind] = sort(ttemp_dose_vec);
% ttemp_count_vec = ttemp_count_vec(ind);
% ttemp_count_sd_vec = ttemp_count_sd_vec(ind);
% x = categorical(string(round(ttemp_dose_vec, 1)));
% x = reordercats(x, string(x));
% 
% figure(403);
% hold on
% h_open = bar(x([3 6 9]), ttemp_count_vec([3 6 9]), 'r', 'LineWidth', 2, ...
%     'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x([3 6 9]), ttemp_count_vec([3 6 9]), ttemp_count_sd_vec([3 6 9]), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% h_dots_peak = bar(x([1 4 7]), ttemp_count_vec([1 4 7]), 'g', 'LineWidth', 2, ...
%     'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x([1 4 7]), ttemp_count_vec([1 4 7]), ttemp_count_sd_vec([1 4 7]), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% h_stripes_peak = bar(x([2 5 8]), ttemp_count_vec([2 5 8]), 'b', ...
%     'LineWidth', 2, 'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x([2 5 8]), ttemp_count_vec([2 5 8]), ttemp_count_sd_vec([2 5 8]), ...
%     'Color', [0 0 0], 'LineStyle', 'none')
% %     lgdstr{k} = [lgdstr{k} '($N=$' num2str() ')']; % Legg til antall kvadranter her for
% hold off
% legend([h_open h_dots_peak h_stripes_peak], lgdstr, ...
%     'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% xlabel('Average quadrat dose (Gy)')
% ylabel('Surviving colony count (per quadrat)')
% ylim([0 Inf])
% grid on
% set(gca, 'FontSize', 18)


% h_vec   = [];
% lgdstr  = {'Peak, Dots', 'Valley, Dots', 'Peak, Stripes', 'Valley, Stripes'};
% color   = {'blue', 'green'};
% facecolor = {[0.78 1.0 0.46], [1.0 0.8 0.8]};
% hatchtype = {'speckle', 'single'};
% temp_dose_vec = [temp_peak_dose_vec temp_valley_dose_vec];
% temp_count_vec = [temp_peak_count_vec temp_valley_count_vec];
% temp_count_sd_vec = [temp_peak_count_sd_vec temp_valley_count_sd_vec];
% x_ctrl = categorical(string(round(temp_dose_vec(1),1)));
% x_ctrl = reordercats(x_ctrl, string(x_ctrl));
% x_open = categorical(string(round(temp_dose_vec(8:10),1)));
% x_open = reordercats(x_open, string(x_open));
% x_peak = categorical(string(round(temp_peak_dose_vec,1)));
% x_peak = reordercats(x_peak, string(x_peak));
% x_valley = categorical(string(round(temp_valley_dose_vec,1)));
% x_valley = reordercats(x_valley, string(x_valley));
% 
% f = figure(401);
% hold on
% h_peak = bar(x_peak((3*k-2):(3*k)), temp_peak_count_vec((3*k-2):(3*k)), ...
%     color{k}, 'LineWidth', 2, 'FaceColor', facecolor{1}, ...
%     'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x_peak((3*k-2):(3*k)), temp_peak_count_vec((3*k-2):(3*k)), ...
%     temp_peak_count_sd_vec((3*k-2):(3*k)), 'Color', [0 0 0], ...
%     'LineStyle', 'none')
% h_valley = bar(x_valley((3*k-2):(3*k)), temp_valley_count_vec((3*k-2):(3*k)), ...
%     color{k}, 'LineWidth', 2, 'FaceColor', facecolor{2}, ...
%     'EdgeColor', [0.2 0.2 0.2], 'FaceAlpha', 0.7);
% errorbar(x_valley((3*k-2):(3*k)), temp_valley_count_vec((3*k-2):(3*k)), ...
%     temp_valley_count_sd_vec((3*k-2):(3*k)), 'Color', [0 0 0], ...
%     'LineStyle', 'none')
% h_vec = [h_vec h_peak h_valley];
% %     lgdstr{k} = [lgdstr{k} '($N=$' num2str() ')']; % Legg til antall
% %     kvadranter her for
% hold off
% legend(h_vec, lgdstr, 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% xlabel('Average quadrat dose (Gy)')
% ylabel('Surviving colony count (per quadrat)')
% ylim([0 Inf])
% grid on
% set(gca, 'FontSize', 18)


% h_vec_1   = [];
% h_vec_2   = [];
% 
% lgdstr  = {'Control', 'Dots', 'Stripes', 'Open'};
% color   = {'black', 'red', 'green', 'blue'};
% 
% for j = 1:length(assay_pattern)
% 
%     N_colonies_vec  = [];
%     dose_avg_vec    = [];
%     SF_vec          = [];
% 
%     for i = 1:length(assay.(assay_pattern{j}).img)
% 
%         dose_avg_vec   = [dose_avg_vec ...
%             assay.(assay_pattern{j}).img{i}.dose_avg];
%         N_colonies_vec = [N_colonies_vec ...
%             assay.(assay_pattern{j}).img{i}.N_colonies];  
% 
%     end
% 
%     SF_vec          = N_colonies_vec ./ N_colonies_ctrl;
% 
%     figure(400);
%     hold on
%     hE = plot(dose_avg_vec, SF_vec, 'o', 'MarkerSize', 6, ...
%         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{j});
%     h_vec_1 = [h_vec_1 hE];
%     hold off
%     xlabel('Average dose (Gy)')
%     ylabel('Surviving fraction')
%     grid on
%     set(gca, 'FontSize', 18)
% 
%     figure(500);
%     hold on
%     hE = plot(dose_avg_vec, N_colonies_vec, 'o', 'MarkerSize', 6, ...
%         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color{j});
%     h_vec_2 = [h_vec_2 hE];
%     hold off
%     xlabel('Average cell flask dose (Gy)')
%     ylabel('Surviving colonies per cell flask')
%     grid on
%     set(gca, 'FontSize', 18)
% 
% 
% end
% legend(h_vec_1, lgdstr, 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% legend(h_vec_2, lgdstr, 'Location', 'NorthEast', 'Interpreter', 'LaTeX')
% xlim([0.375 5.8])
% ylim([0 1.2])

%% K-fold cross-validation for Poisson regression of quadrat data (training)

close all

path        = 'C:\Users\delmo\Desktop\New folder1\';

% Generate K-fold indices
K = 4;
cv_indices = crossvalind('Kfold', 12, K);

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
eta_LQ   = []; eta_MLQ_R   = []; eta_MLQ_A   = []; eta_MLQ_R2   = []; eta_MLQ_A2   = [];

for k = 1:K

    test_indices  = (cv_indices == k);
    train_indices = ~test_indices;
    test_indices  = find(test_indices);
    train_indices = find(train_indices);

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

            assay_test.(assay_pattern{j})  = assay.(assay_pattern{j});
            assay_train.(assay_pattern{j}) = assay.(assay_pattern{j});

        end

    end

    % Design variable matrix X
    X_D     = [];
    X_R     = [];
    X_A     = [];
    Y_CC    = [];
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
        X_D     = [X_D      dose_vec.(assay_pattern{j})];
        X_R     = [X_R      peakdist_vec.(assay_pattern{j})];
        X_G     = [X_G      grad_vec.(assay_pattern{j})];
        Y_CC    = [Y_CC     count_vec.(assay_pattern{j})];
        if any(strcmp(assay_pattern{j}, {'Control', 'Open'}))
            X_Open  = [X_Open   dose_vec.(assay_pattern{j})];
            Y_Open  = [Y_Open   count_vec.(assay_pattern{j})];
        end

    end

%     [X_D, index_sort]       = sort(X_D);
%     X_R                     = X_R(index_sort);
%     X_A                     = X_A(index_sort);
%     Y_CC                    = Y_CC(index_sort);
%     [X_Open, index_sort]    = sort(X_Open);
%     Y_Open                  = Y_Open(index_sort);

    % Perform Poisson regression
    % alpha = 0.24 +- 0.02 Gy-1, beta = 0.019 +- 0.002 Gy-2
    tbl = table(X_D', X_R', X_G', Y_CC', 'VariableNames', {'D', 'R', 'G', 'Y_CC'});
    tbl_Open = table(X_Open', Y_Open', 'VariableNames', {'D_Open', 'Y_CC_Open'});

    mdl_LQ = fitglm(tbl, 'Y_CC ~ 1 + D + D^2', 'Link', 'log', ...
        'Distribution', 'poisson', 'Intercept', true);
    mdl_MLQ_R = fitglm(tbl, 'Y_CC ~ 1 + D + D^2 + R', 'Link', 'log', ...
        'Distribution', 'poisson', 'Intercept', true);
    mdl_MLQ_G = fitglm(tbl, 'Y_CC ~ 1 + D + D^2 + G', 'Link', 'log', ...
        'Distribution', 'poisson', 'Intercept', true);
    mdl_MLQ_DG = fitglm(tbl, 'Y_CC ~ 1 + D + D^2 + D:G', 'Link', 'log', ...
        'Distribution', 'poisson', 'Intercept', true);
    mdl_Open = fitglm(tbl_Open, 'Y_CC_Open ~ 1 + D_Open + D_Open^2', ...
        'Link', 'log', 'Distribution', 'poisson', 'Intercept', true);
    %     [mdl_Open.ModelCriterion.AIC mdl_Open.ModelCriterion.BIC; ...
    %         mdl_LQ.ModelCriterion.AIC mdl_LQ.ModelCriterion.BIC; ...
    %         mdl_MLQ_R.ModelCriterion.AIC mdl_MLQ_R.ModelCriterion.BIC]
    [mdl_LQ.ModelCriterion.AIC mdl_LQ.ModelCriterion.BIC; ...
        mdl_MLQ_R.ModelCriterion.AIC mdl_MLQ_R.ModelCriterion.BIC; ...
        mdl_MLQ_G.ModelCriterion.AIC mdl_MLQ_G.ModelCriterion.BIC; ...
        mdl_MLQ_DG.ModelCriterion.AIC mdl_MLQ_DG.ModelCriterion.BIC]
    %mdl_MLQ_DR.ModelCriterion.AIC mdl_MLQ_DR.ModelCriterion.BIC; ...
    %mdl_MLQ_DA.ModelCriterion.AIC mdl_MLQ_DA.ModelCriterion.BIC]

%     G = mdl_LQ.Deviance - mdl_MLQ_R.Deviance;
%     p = 1 - chi2cdf(G, 1)

    % Create data points for prediction
    D_fit            = linspace(min(X_D(:)), max(X_D(:)), 101);
    R_fit            = linspace(min(X_R(:)), max(X_R(:)), 101);
    G_fit            = linspace(min(X_G(:)), max(X_G(:)), 101);
    [~,      R_mesh] = meshgrid(D_fit, R_fit);
    [D_mesh, G_mesh] = meshgrid(D_fit, G_fit);

    % LQ prediction
    lambda0_LQ(k) = mdl_LQ.Coefficients.Estimate(1);
    alpha_LQ(k)   = mdl_LQ.Coefficients.Estimate(2);
    beta_LQ(k)    = mdl_LQ.Coefficients.Estimate(3);
    eta_LQ(k)     = 0.0;
    Y_CC_fit_LQ   = exp(lambda0_LQ(k) + alpha_LQ(k).*D_mesh + ...
        beta_LQ(k).*(D_mesh.^2) + eta_LQ(k).*R_mesh);

    % MLQ (R) prediciton
    lambda0_MLQ_R(k) = mdl_MLQ_R.Coefficients.Estimate(1);
    alpha_MLQ_R(k)   = mdl_MLQ_R.Coefficients.Estimate(2);
    beta_MLQ_R(k)    = mdl_MLQ_R.Coefficients.Estimate(4);
    eta_MLQ_R(k)     = mdl_MLQ_R.Coefficients.Estimate(3);
    Y_CC_fit_MLQ_R = exp(lambda0_MLQ_R(k) + alpha_MLQ_R(k).*D_mesh + ...
        beta_MLQ_R(k).*(D_mesh.^2) + eta_MLQ_R(k).*(R_mesh));

    % MLQ (A) prediciton
    lambda0_MLQ_G(k) = mdl_MLQ_G.Coefficients.Estimate(1);
    alpha_MLQ_G(k)   = mdl_MLQ_G.Coefficients.Estimate(2);
    beta_MLQ_G(k)    = mdl_MLQ_G.Coefficients.Estimate(4);
    eta_MLQ_G(k)     = mdl_MLQ_G.Coefficients.Estimate(3);
    Y_CC_fit_MLQ_G = exp(lambda0_MLQ_G(k) + alpha_MLQ_G(k).*D_mesh + ...
        beta_MLQ_G(k).*(D_mesh.^2) + eta_MLQ_G(k).*(G_mesh));

    % MLQ (R^2) prediciton
    lambda0_MLQ_DG(k) = mdl_MLQ_DG.Coefficients.Estimate(1);
    alpha_MLQ_DG(k)   = mdl_MLQ_DG.Coefficients.Estimate(2);
    beta_MLQ_DG(k)    = mdl_MLQ_DG.Coefficients.Estimate(4);
    eta_MLQ_DG(k)     = mdl_MLQ_DG.Coefficients.Estimate(3);
    Y_CC_fit_MLQ_DG = exp(lambda0_MLQ_DG(k) + alpha_MLQ_DG(k).*D_mesh + ...
        beta_MLQ_DG(k).*(D_mesh.^2) + eta_MLQ_DG(k).*(G_mesh));

    % Plot prediction surface
    plot3Dsurf(X_D, X_R, Y_CC, D_mesh, R_mesh, Y_CC_fit_LQ, Y_CC_fit_MLQ_R, ...
        'Distance to peak (cm)', ...
        'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\deltaR)})', ...
        fullfile(path, sprintf('SurfPlot_LQvsMLQ_R_%s.png', num2str(k))))
    plot3Dsurf(X_D, X_G, Y_CC, D_mesh, G_mesh, Y_CC_fit_LQ, Y_CC_fit_MLQ_G, ...
        'Irradiation fraction', ...
        'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\delta\bar{\dot{D}})})', ...
        fullfile(path, sprintf('SurfPlot_LQvsMLQ_G_%s.png', num2str(k))))
    plot3Dsurf(X_D, X_R, Y_CC, D_mesh, R_mesh, Y_CC_fit_LQ, Y_CC_fit_MLQ_DG, ...
        'Distance to peak (cm)', ...
        'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\deltaD\bar{\dot{D}})})', ...
        fullfile(path, sprintf('SurfPlot_LQvsMLQ_DG_%s.png', num2str(k))))

    % Test
    h_vec                       = [];
    dose_quadrat_pred           = [];
    count_quadrat_LQ            = [];
    count_quadrat_MLQ_R         = [];
    count_quadrat_MLQ_A         = [];
    count_quadrat_MLQ_R2        = [];
    count_peak_quadrat_LQ       = [];
    count_peak_quadrat_MLQ_R    = [];
    count_peak_quadrat_MLQ_G    = [];
    count_peak_quadrat_MLQ_DG   = [];
    count_valley_quadrat_LQ     = [];
    count_valley_quadrat_MLQ_R  = [];
    count_valley_quadrat_MLQ_G  = [];
    count_valley_quadrat_MLQ_DG = [];

    for j = 1:length(assay_pattern)

        dose_quadrat_obs            = [];
        dose_peak_quadrat_obs       = []; 
        dose_valley_quadrat_obs     = [];
        count_quadrat_obs           = []; 
        count_peak_quadrat_obs      = [];
        count_valley_quadrat_obs    = [];

        for i = 1:length(assay_test.(assay_pattern{j}).img)

            str_dose = num2str(assay_test.(assay_pattern{j}).img{i}.dose);

            %img_dose_quadrat       = assay_test.(assay_pattern{j}).img{i}.img_dose_quadrat;
            %img_count_quadrat      = assay_test.(assay_pattern{j}).img{i}.img_count_quadrat;
            temp_img_dose          = assay_test.(assay_pattern{j}).img{i}.img_dose;
            temp_img_peakdist      = assay_test.(assay_pattern{j}).img{i}.img_peakdist;
            temp_img_grad          = assay_test.(assay_pattern{j}).img{i}.img_grad;
            bw_quadrat             = assay_test.(assay_pattern{j}).img{i}.bw_quadrat;
            if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
                bw_peak_quadrat   = assay.(assay_pattern{j}).img{i}.bw_peak_quadrat;
                bw_valley_quadrat = assay.(assay_pattern{j}).img{i}.bw_valley_quadrat;
            end

            % Observation
            PV_dose = assay_test.(assay_pattern{j}).img{i}.PV.dose_quadrat;
            PV_count = assay_test.(assay_pattern{j}).img{i}.PV.count_quadrat;
            dose_quadrat_obs = [dose_quadrat_obs mean(PV_dose)];
            dose_quadrat_pred = [dose_quadrat_pred mean(PV_dose)];
            count_quadrat_obs = [count_quadrat_obs mean(PV_count)];
            if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
                dose_peak_quadrat_obs = [dose_peak_quadrat_obs ...
                    mean(assay_test.(assay_pattern{j}).img{i}.PV.dose_peak_quadrat)];
                count_peak_quadrat_obs = [count_peak_quadrat_obs ...
                    mean(assay_test.(assay_pattern{j}).img{i}.PV.count_peak_quadrat)];
                dose_valley_quadrat_obs = [dose_valley_quadrat_obs ...
                    mean(assay_test.(assay_pattern{j}).img{i}.PV.dose_valley_quadrat)];
                count_valley_quadrat_obs = [count_valley_quadrat_obs ...
                    mean(assay_test.(assay_pattern{j}).img{i}.PV.count_valley_quadrat)];
            end

            % Prediction
            img_count_LQ    = LQ(temp_img_dose, lambda0_LQ(k), alpha_LQ(k), beta_LQ(k));
            img_count_MLQ_R = MLQ(temp_img_dose, temp_img_peakdist, ...
                lambda0_MLQ_R(k), alpha_MLQ_R(k), beta_MLQ_R(k), eta_MLQ_R(k), 'lin');
            img_count_MLQ_G = MLQ(temp_img_dose, temp_img_grad, ...
                lambda0_MLQ_G(k), alpha_MLQ_G(k), beta_MLQ_G(k), eta_MLQ_G(k), 'lin');
            img_count_MLQ_DG = MLQ(temp_img_dose, temp_img_peakdist, ...
                lambda0_MLQ_DG(k), alpha_MLQ_DG(k), beta_MLQ_DG(k), eta_MLQ_DG(k), 'sq');

            % Get all relevant quadrat colony count data
            PV_count_LQ    = statsModel(img_count_LQ, bw_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
            PV_count_MLQ_R = statsModel(img_count_MLQ_R, bw_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
            PV_count_MLQ_G = statsModel(img_count_MLQ_G, bw_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
            PV_count_MLQ_DG = statsModel(img_count_MLQ_DG, bw_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
            count_quadrat_LQ     = [count_quadrat_LQ mean(PV_count_LQ)];
            count_quadrat_MLQ_R  = [count_quadrat_MLQ_R mean(PV_count_MLQ_R)];
            count_quadrat_MLQ_G  = [count_quadrat_MLQ_G mean(PV_count_MLQ_G)];
            count_quadrat_MLQ_DG = [count_quadrat_MLQ_DG mean(PV_count_MLQ_DG)];
            if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
                PV_count_LQ_peak       = statsModel(img_count_LQ, bw_peak_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_MLQ_R_peak    = statsModel(img_count_MLQ_R, bw_peak_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_MLQ_G_peak    = statsModel(img_count_MLQ_G, bw_peak_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_MLQ_DG_peak   = statsModel(img_count_MLQ_DG, bw_peak_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_LQ_valley     = statsModel(img_count_LQ, bw_valley_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_MLQ_R_valley  = statsModel(img_count_MLQ_R, bw_valley_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_MLQ_G_valley  = statsModel(img_count_MLQ_G, bw_valley_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                PV_count_MLQ_DG_valley = statsModel(img_count_MLQ_DG, bw_valley_quadrat, n_quadrat_x, n_quadrat_y, dxdy_px);
                count_peak_quadrat_LQ       = [count_peak_quadrat_LQ mean(PV_count_LQ_peak)];
                count_peak_quadrat_MLQ_R    = [count_peak_quadrat_MLQ_R mean(PV_count_MLQ_R_peak)];
                count_peak_quadrat_MLQ_G    = [count_peak_quadrat_MLQ_A mean(PV_count_MLQ_G_peak)];
                count_peak_quadrat_MLQ_DG   = [count_peak_quadrat_MLQ_R2 mean(PV_count_MLQ_DG_peak)];
                count_valley_quadrat_LQ     = [count_valley_quadrat_LQ mean(PV_count_LQ_valley)];
                count_valley_quadrat_MLQ_R  = [count_valley_quadrat_MLQ_R mean(PV_count_MLQ_R_valley)];
                count_valley_quadrat_MLQ_G  = [count_valley_quadrat_MLQ_G mean(PV_count_MLQ_G_valley)];
                count_valley_quadrat_MLQ_DG = [count_valley_quadrat_MLQ_DG mean(PV_count_MLQ_DG_valley)];
            end

            % Store results
            if isfield(assay_test.(assay_pattern{j}), ['Dose' str_dose 'Gy'])
                sz = size(assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_LQ, 3);
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_dose(:,:,sz+1)            = temp_img_dose;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).bw_quadrat(:,:,sz+1)          = bw_quadrat;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_LQ(:,:,sz+1)        = img_count_LQ;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_MLQ_R(:,:,sz+1)     = img_count_MLQ_R;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_MLQ_G(:,:,sz+1)     = img_count_MLQ_G;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_MLQ_DG(:,:,sz+1)    = img_count_MLQ_DG;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.dose_quadrat{sz+1}         = PV_dose;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_quadrat{sz+1}        = PV_count;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_LQ_quadrat{sz+1}     = PV_count_LQ;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_MLQ_R_quadrat{sz+1}  = PV_count_MLQ_R;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_MLQ_G_quadrat{sz+1}  = PV_count_MLQ_G;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_MLQ_DG_quadrat{sz+1} = PV_count_MLQ_DG;
            else
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_dose(:,:,1)            = temp_img_dose;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).bw_quadrat(:,:,1)          = bw_quadrat;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_LQ(:,:,1)        = img_count_LQ;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_MLQ_R(:,:,1)     = img_count_MLQ_R;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_MLQ_G(:,:,1)     = img_count_MLQ_G;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).img_count_MLQ_DG(:,:,1)    = img_count_MLQ_DG;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.dose_quadrat{1}         = PV_dose;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_quadrat{1}        = PV_count;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_LQ_quadrat{1}     = PV_count_LQ;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_MLQ_R_quadrat{1}  = PV_count_MLQ_R;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_MLQ_G_quadrat{1}  = PV_count_MLQ_G;
                assay_test.(assay_pattern{j}).(['Dose' str_dose 'Gy']).PV.count_MLQ_DG_quadrat{1} = PV_count_MLQ_DG;
            end

        end

%         figure(400*k);
%         hold on
%         h = plot(dose_quadrat_obs, count_quadrat_obs, 'o', ...
%             'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
%             'MarkerFaceColor', color{j});
%         h_vec = [h_vec h];
%         if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
%             h_peak = plot(dose_peak_quadrat_obs, count_peak_quadrat_obs, ...
%                 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
%                 'MarkerFaceColor', color{j});
%             h_valley = plot(dose_valley_quadrat_obs, count_valley_quadrat_obs, ...
%                 '^', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
%                 'MarkerFaceColor', color{j});
%             h_vec = [h_vec h_peak h_valley];
%         end

    end

%     Y_CC_fit_LQ     = exp(lambda0_LQ + alpha_LQ.*D_fit + beta_LQ.*(D_fit.^2));
%     Y_CC_fit_MLQ_R  = [];
%     Y_CC_fit_MLQ_A  = [];
%     Y_CC_fit_MLQ_R2 = [];
%     Y_CC_fit_MLQ_A2 = [];
%     for ii = 1:length(R_fit)
%         temp_Y_CC_fit_MLQ_R = exp(lambda0_MLQ_R + alpha_MLQ_R.*D_fit + ...
%             beta_MLQ_R.*(D_fit.^2) + eta_MLQ_R*(R_fit(ii)));
%         temp_Y_CC_fit_MLQ_A = exp(lambda0_MLQ_A + alpha_MLQ_A.*D_fit + ...
%             beta_MLQ_A.*(D_fit.^2) + eta_MLQ_A.*(A_fit(ii)));
%         temp_Y_CC_fit_MLQ_R2 = exp(lambda0_MLQ_R2 + alpha_MLQ_R2.*D_fit + ...
%             beta_MLQ_R2.*(D_fit.^2) + eta_MLQ_R2.*(R_fit(ii).^2));
%         temp_Y_CC_fit_MLQ_A2 = exp(lambda0_MLQ_A2 + alpha_MLQ_A2.*D_fit + ...
%             beta_MLQ_A2.*(D_fit.^2) + eta_MLQ_A2.*(A_fit(ii).^2));
%         Y_CC_fit_MLQ_R  = [Y_CC_fit_MLQ_R;  temp_Y_CC_fit_MLQ_R];
%         Y_CC_fit_MLQ_A  = [Y_CC_fit_MLQ_A;  temp_Y_CC_fit_MLQ_A];
%         Y_CC_fit_MLQ_R2 = [Y_CC_fit_MLQ_R2; temp_Y_CC_fit_MLQ_R2];
%         Y_CC_fit_MLQ_A2 = [Y_CC_fit_MLQ_A2; temp_Y_CC_fit_MLQ_A2];
%     end
%     profile_MLQ_R  = estimateCIBand(Y_CC_fit_MLQ_R, D_fit, significance);
%     profile_MLQ_A  = estimateCIBand(Y_CC_fit_MLQ_A, D_fit, significance);
%     profile_MLQ_R2 = estimateCIBand(Y_CC_fit_MLQ_R2, D_fit, significance);
%     profile_MLQ_A2 = estimateCIBand(Y_CC_fit_MLQ_A2, D_fit, significance);
% 
%     h_LQ = plot(D_fit, Y_CC_fit_LQ, linespec{1});
%     h_95CImean = plot(profile_MLQ_R.x, profile_MLQ_R.y_avg, ...
%         linespec{2}, 'LineWidth', 1.0);
%     h_95CIband = fill(profile_MLQ_R.CI95_x, profile_MLQ_R.CI95_value, ...
%         color_CI95{2}, 'FaceColor', facecolor_CI95{2}, 'EdgeColor', 'none', ...
%         'FaceAlpha', 0.5);
%     hold off
%     h_vec = [h_vec h_LQ h_95CImean h_95CIband];
%     legend(h_vec, lgdstr, 'Location', 'NorthEast')
%     xlim([0 10])
%     xlabel('Average dose (Gy)')
%     ylabel('Surviving colony count (per quadrat)')
%     grid on
%     set(gca, 'FontSize', 18)

end



%% Plot

close all

PV_struct = {};
step = [0 0.1 0.1 0.001];

for j = 2:length(assay_pattern)

    PV_struct.(assay_pattern{j}).dose_array         = [];
    PV_struct.(assay_pattern{j}).LQ.SF.avg          = [];
    PV_struct.(assay_pattern{j}).LQ.SF.SD           = [];
    PV_struct.(assay_pattern{j}).MLQ_R.SF.avg       = [];
    PV_struct.(assay_pattern{j}).MLQ_R.SF.SD        = [];
    PV_struct.(assay_pattern{j}).MLQ_A.SF.avg       = [];
    PV_struct.(assay_pattern{j}).MLQ_A.SF.SD        = [];
    PV_struct.(assay_pattern{j}).MLQ_R2.SF.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_R2.SF.SD       = [];
    PV_struct.(assay_pattern{j}).MLQ_A2.SF.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_A2.SF.SD       = [];
    PV_struct.(assay_pattern{j}).MLQ_DR.SF.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_DR.SF.SD       = [];
    PV_struct.(assay_pattern{j}).MLQ_DA.SF.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_DA.SF.SD       = [];
    PV_struct.(assay_pattern{j}).LQ.MAE.avg         = [];
    PV_struct.(assay_pattern{j}).MLQ_R.MAE.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_A.MAE.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_R2.MAE.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_A2.MAE.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_DR.MAE.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_DA.MAE.avg     = [];
    PV_struct.(assay_pattern{j}).LQ.RMSE.avg        = [];
    PV_struct.(assay_pattern{j}).MLQ_R.RMSE.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_A.RMSE.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_R2.RMSE.avg    = [];
    PV_struct.(assay_pattern{j}).MLQ_A2.RMSE.avg    = [];
    PV_struct.(assay_pattern{j}).MLQ_DR.RMSE.avg    = [];
    PV_struct.(assay_pattern{j}).MLQ_DA.RMSE.avg    = [];
    PV_struct.(assay_pattern{j}).LQ.LogLik.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_R.LogLik.avg   = [];
    PV_struct.(assay_pattern{j}).MLQ_A.LogLik.avg   = [];
    PV_struct.(assay_pattern{j}).MLQ_R2.LogLik.avg  = [];
    PV_struct.(assay_pattern{j}).MLQ_A2.LogLik.avg  = [];
    PV_struct.(assay_pattern{j}).MLQ_DR.LogLik.avg  = [];
    PV_struct.(assay_pattern{j}).MLQ_DA.LogLik.avg  = [];
    PV_struct.(assay_pattern{j}).LQ.Rsq.avg         = [];
    PV_struct.(assay_pattern{j}).MLQ_R.Rsq.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_A.Rsq.avg      = [];
    PV_struct.(assay_pattern{j}).MLQ_R2.Rsq.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_A2.Rsq.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_DR.Rsq.avg     = [];
    PV_struct.(assay_pattern{j}).MLQ_DA.Rsq.avg     = [];

    for i = 1:length(assay_dose)

%         temp_img_SF_LQ      = ...
%             sum(assay_test.(assay_pattern{j}).(assay_dose{i}).img_SF_LQ, 3) ./ ...
%             size(assay_test.(assay_pattern{j}).(assay_dose{i}).img_SF_LQ, 3);
%         temp_img_SF_MLQ     = ...
%             sum(assay_test.(assay_pattern{j}).(assay_dose{i}).img_SF_MLQ, 3) ./ ...
%             size(assay_test.(assay_pattern{j}).(assay_dose{i}).img_SF_MLQ, 3);
%         temp_img_RPD_LQ     = ...
%             sum(assay_test.(assay_pattern{j}).(assay_dose{i}).img_diff_LQ_quadrat, 3) ./ ...
%             size(assay_test.(assay_pattern{j}).(assay_dose{i}).img_diff_LQ_quadrat, 3);
%         temp_img_RPD_MLQ    = ...
%             sum(assay_test.(assay_pattern{j}).(assay_dose{i}).img_diff_MLQ_quadrat, 3) ./ ...
%             size(assay_test.(assay_pattern{j}).(assay_dose{i}).img_diff_MLQ_quadrat, 3);

        temp_PV_dose      = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.dose_quadrat{:});
        temp_PV_SF        = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_quadrat{:});
        temp_PV_SF_LQ     = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_LQ_quadrat{:});
        temp_PV_SF_MLQ_R  = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_MLQ_R_quadrat{:});
        temp_PV_SF_MLQ_A  = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_MLQ_A_quadrat{:});
        temp_PV_SF_MLQ_R2 = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_MLQ_R2_quadrat{:});
        temp_PV_SF_MLQ_A2 = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_MLQ_A2_quadrat{:});
%         temp_PV_SF_MLQ_DR = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_MLQ_DR_quadrat{:});
%         temp_PV_SF_MLQ_DA = vertcat(assay_test.(assay_pattern{j}).(assay_dose{i}).PV.SF_MLQ_DA_quadrat{:});

        %         [temp_PV_dose, index_sort] = sort(temp_PV_dose);
        %         temp_PV_SF_LQ              = temp_PV_SF_LQ(index_sort);
        %         temp_PV_SF_MLQ             = temp_PV_SF_MLQ(index_sort);

        % Compute mean absolute error (MAE), root-mean-squared error (RMSE)
        % and R-squared performance
        % PV_dose.(assay_pattern{j}).dose(i)             = temp_PV_dose;
        PV_struct.(assay_pattern{j}).LQ.SF.avg(i)       = mean(temp_PV_SF_LQ);
        PV_struct.(assay_pattern{j}).LQ.SF.SD(i)        = std(temp_PV_SF_LQ);
        PV_struct.(assay_pattern{j}).MLQ_R.SF.avg(i)    = mean(temp_PV_SF_MLQ_R);
        PV_struct.(assay_pattern{j}).MLQ_R.SF.SD(i)     = std(temp_PV_SF_MLQ_R);
        PV_struct.(assay_pattern{j}).MLQ_A.SF.avg(i)    = mean(temp_PV_SF_MLQ_A);
        PV_struct.(assay_pattern{j}).MLQ_A.SF.SD(i)     = std(temp_PV_SF_MLQ_A);
        PV_struct.(assay_pattern{j}).MLQ_R2.SF.avg(i)   = mean(temp_PV_SF_MLQ_R2);
        PV_struct.(assay_pattern{j}).MLQ_R2.SF.SD(i)    = std(temp_PV_SF_MLQ_R2);
        PV_struct.(assay_pattern{j}).MLQ_A2.SF.avg(i)   = mean(temp_PV_SF_MLQ_A2);
        PV_struct.(assay_pattern{j}).MLQ_A2.SF.SD(i)    = std(temp_PV_SF_MLQ_A2);
        PV_struct.(assay_pattern{j}).MLQ_DR.SF.avg(i)   = mean(temp_PV_SF_MLQ_DR);
        PV_struct.(assay_pattern{j}).MLQ_DR.SF.SD(i)    = std(temp_PV_SF_MLQ_DR);
        PV_struct.(assay_pattern{j}).MLQ_DA.SF.avg(i)   = mean(temp_PV_SF_MLQ_DA);
        PV_struct.(assay_pattern{j}).MLQ_DA.SF.SD(i)    = std(temp_PV_SF_MLQ_DA);
        PV_struct.(assay_pattern{j}).LQ.MAE.avg(i)      = mean(abs(temp_PV_SF_LQ - temp_PV_SF));
        PV_struct.(assay_pattern{j}).MLQ_R.MAE.avg(i)   = mean(abs(temp_PV_SF_MLQ_R - temp_PV_SF));
        PV_struct.(assay_pattern{j}).MLQ_A.MAE.avg(i)   = mean(abs(temp_PV_SF_MLQ_A - temp_PV_SF));
        PV_struct.(assay_pattern{j}).MLQ_R2.MAE.avg(i)  = mean(abs(temp_PV_SF_MLQ_R2 - temp_PV_SF));
        PV_struct.(assay_pattern{j}).MLQ_A2.MAE.avg(i)  = mean(abs(temp_PV_SF_MLQ_A2 - temp_PV_SF));
        PV_struct.(assay_pattern{j}).MLQ_DR.MAE.avg(i)  = mean(abs(temp_PV_SF_MLQ_DR - temp_PV_SF));
        PV_struct.(assay_pattern{j}).MLQ_DA.MAE.avg(i)  = mean(abs(temp_PV_SF_MLQ_DA - temp_PV_SF));
        PV_struct.(assay_pattern{j}).LQ.RMSE.avg(i)     = sqrt(mean((temp_PV_SF_LQ - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).MLQ_R.RMSE.avg(i)  = sqrt(mean((temp_PV_SF_MLQ_R - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).MLQ_A.RMSE.avg(i)  = sqrt(mean((temp_PV_SF_MLQ_A - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).MLQ_R2.RMSE.avg(i) = sqrt(mean((temp_PV_SF_MLQ_R2 - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).MLQ_A2.RMSE.avg(i) = sqrt(mean((temp_PV_SF_MLQ_A2 - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).MLQ_DR.RMSE.avg(i) = sqrt(mean((temp_PV_SF_MLQ_DR - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).MLQ_DA.RMSE.avg(i) = sqrt(mean((temp_PV_SF_MLQ_DA - temp_PV_SF).^2));
        PV_struct.(assay_pattern{j}).LQ.LogLik.avg(i)     = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_LQ);
        PV_struct.(assay_pattern{j}).MLQ_R.LogLik.avg(i)  = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_MLQ_R);
        PV_struct.(assay_pattern{j}).MLQ_A.LogLik.avg(i)  = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_MLQ_A);
        PV_struct.(assay_pattern{j}).MLQ_R2.LogLik.avg(i) = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_MLQ_R2);
        PV_struct.(assay_pattern{j}).MLQ_A2.LogLik.avg(i) = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_MLQ_A2);
        PV_struct.(assay_pattern{j}).MLQ_DR.LogLik.avg(i) = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_MLQ_DR);
        PV_struct.(assay_pattern{j}).MLQ_DA.LogLik.avg(i) = estimateLogLikelihoodPoisson(temp_PV_SF, temp_PV_SF_MLQ_DA);
        PV_struct.(assay_pattern{j}).LQ.Rsq.avg(i)      = corr(temp_PV_SF_LQ, temp_PV_SF)^2;
        PV_struct.(assay_pattern{j}).MLQ_R.Rsq.avg(i)   = corr(temp_PV_SF_MLQ_R, temp_PV_SF)^2;
        PV_struct.(assay_pattern{j}).MLQ_A.Rsq.avg(i)   = corr(temp_PV_SF_MLQ_A, temp_PV_SF)^2;
        PV_struct.(assay_pattern{j}).MLQ_R2.Rsq.avg(i)  = corr(temp_PV_SF_MLQ_R2, temp_PV_SF)^2;
        PV_struct.(assay_pattern{j}).MLQ_A2.Rsq.avg(i)  = corr(temp_PV_SF_MLQ_A2, temp_PV_SF)^2;
        PV_struct.(assay_pattern{j}).MLQ_DR.Rsq.avg(i)  = corr(temp_PV_SF_MLQ_DR, temp_PV_SF)^2;
        PV_struct.(assay_pattern{j}).MLQ_DA.Rsq.avg(i)  = corr(temp_PV_SF_MLQ_DA, temp_PV_SF)^2;

        %         h = figure();
        %         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %         hold on
        %         histogram2(temp_PV_dose, temp_PV_SF_LQ, 'NumBins', ...
        %             [60 round((max(temp_PV_SF_LQ)-min(temp_PV_SF_LQ))/step(j))], ...
        %             'FaceColor', 'green', 'EdgeColor', 'k', 'FaceAlpha', 0.5)
        %         histogram2(temp_PV_dose, temp_PV_SF_MLQ, 'NumBins', ...
        %             [60 round((max(temp_PV_SF_MLQ)-min(temp_PV_SF_MLQ))/step(j))], ...
        %             'FaceColor', 'blue', 'EdgeColor', 'k', 'FaceAlpha', 0.5)
        %         xlim([min(temp_PV_dose) max(temp_PV_dose)])
        %         ylim([min(min(temp_PV_SF_LQ, temp_PV_SF_MLQ)) ...
        %             max(max(temp_PV_SF_LQ, temp_PV_SF_MLQ))])
        %         if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
        %             line([peak_dose.(assay_pattern{j}).(assay_dose{i}), ...
        %                 peak_dose.(assay_pattern{j}).(assay_dose{i})], ylim, [0,0], ...
        %                 'LineStyle', '--', 'LineWidth', 2.0, 'Color', 'red')
        %             line([valley_dose.(assay_pattern{j}).(assay_dose{i}), ...
        %                 valley_dose.(assay_pattern{j}).(assay_dose{i})], ylim, [0,0], ...
        %                 'LineStyle', '--', 'LineWidth', 2.0, 'Color', 'cyan')
        %         end
        %         hold off
        %         xlabel('Dose (Gy)');
        %         ylabel('Predicted colony count per quadrat');
        %         zlabel('Number of quadrats')
        %         grid on; box on; camlight; view(40,35)
        %         if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
        %             lgd = legend('LQ', 'MLQ', 'Peak dose', 'Valley dose', ...
        %                 'NumColumns', 2, 'Location', 'NorthEast');
        %         else
        %             lgd = legend('LQ', 'MLQ', 'Location', 'NorthEast');
        %         end
        %         title(lgd, [pattern{j} ' irradiaton; ' dose_str{i} ' nominal dose'])
        %         set(gca, 'FontSize', 16)
        %
        %         h = figure();
        %         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %         hold on
        %         histogram2(temp_PV_dose, temp_PV_RPD_LQ, 'NumBins', ...
        %             [60 round((max(temp_PV_RPD_LQ)-min(temp_PV_RPD_LQ))/step(j))], ...
        %             'FaceColor', 'green', 'EdgeColor', 'k', 'FaceAlpha', 0.5)
        %         histogram2(temp_PV_dose, temp_PV_RPD_MLQ, 'NumBins', ...
        %             [60 round((max(temp_PV_RPD_MLQ)-min(temp_PV_RPD_MLQ))/step(j))], ...
        %             'FaceColor', 'blue', 'EdgeColor', 'k', 'FaceAlpha', 0.5)
        %         xlim([min(temp_PV_dose) max(temp_PV_dose)])
        %         ylim([min(min(temp_PV_RPD_LQ, temp_PV_RPD_MLQ)) ...
        %             max(max(temp_PV_RPD_LQ, temp_PV_RPD_MLQ))])
        %         if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
        %             line([peak_dose.(assay_pattern{j}).(assay_dose{i}), ...
        %                 peak_dose.(assay_pattern{j}).(assay_dose{i})], ylim, [0,0], ...
        %                 'LineStyle', '--', 'LineWidth', 2.0, 'Color', 'red')
        %             line([valley_dose.(assay_pattern{j}).(assay_dose{i}), ...
        %                 valley_dose.(assay_pattern{j}).(assay_dose{i})], ylim, [0,0], ...
        %                 'LineStyle', '--', 'LineWidth', 2.0, 'Color', 'cyan')
        %         end
        %         hold off
        %         xlabel('Dose (Gy)');
        %         ylabel('Relative deviation');
        %         zlabel('Number of quadrats')
        %         grid on; box on; camlight; view(40,35)
        %         if any(strcmp(assay_pattern{j}, {'GRID_Stripes', 'GRID_Dots'}))
        %             lgd = legend('LQ', 'MLQ', 'Peak dose', 'Valley dose', ...
        %                 'NumColumns', 2, 'Location', 'NorthEast');
        %         else
        %             lgd = legend('LQ', 'MLQ', 'Location', 'NorthEast');
        %         end
        %         title(lgd, [pattern{j} ' irradiaton; ' dose_str{i} ' nominal dose'])
        %         set(gca, 'FontSize', 16)

        %         h = figure();
        %         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %         imshow(temp_img_SF_LQ, [0 max([ceil(max(max(temp_img_SF_LQ))) ...
        %             ceil(max(max(temp_img_SF_MLQ)))])], 'colormap', jet(4096))
        %         shading interp
        %         c = colorbar;
        %         c.Label.String = 'Colony count';
        %
        %         h = figure();
        %         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %         imshow(temp_img_SF_MLQ, [0 max([ceil(max(max(temp_img_SF_LQ))) ...
        %             ceil(max(max(temp_img_SF_MLQ)))])], 'colormap', jet(4096))
        %         shading interp
        %         c = colorbar;
        %         c.Label.String = 'Colony count';
        %
        %         plot2Dmap(temp_img_RPD_LQ, 47.*px_size_assay, ...
        %             [min(min(temp_img_RPD_LQ(:)), min(temp_img_RPD_MLQ(:))) ...
        %             0.95*max(max(temp_img_RPD_LQ(:)), max(temp_img_RPD_MLQ(:)))], ...
        %             'Relative percentage deviation (%)', ...
        %             ['LQ; ' dose_str{i} ' nominal dose'])
        %
        %         plot2Dmap(temp_img_RPD_MLQ, 47.*px_size_assay, ...
        %             [min(min(temp_img_RPD_LQ(:)), min(temp_img_RPD_MLQ(:))) ...
        %             0.95*max(max(temp_img_RPD_LQ(:)), max(temp_img_RPD_MLQ(:)))], ...
        %             'Relative percentage deviation (%)', ...
        %             ['MLQ; ' dose_str{i} ' nominal dose'])

    end

    h_pred = figure();
    set(h_pred, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    hE_LQ = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-og', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');
    hold on
    hE_MLQ_R = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-ob', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
    hE_MLQ_A = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-or', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
    hE_MLQ_R2 = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-oc', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c');
    hE_MLQ_A2 = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-om', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm');
    hE_MLQ_DR = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-oy', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y');
    hE_MLQ_DA = errorbar(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.SF.avg, ...
        PV_struct.(assay_pattern{j}).LQ.SF.SD, '-ok', 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    hold off
    xlabel('Dose (Gy)');
    ylabel('Prediction (colony count per quadrat)');
    xlim([min(unique(dose_assay))-1 max(unique(dose_assay))+1])
%     ylim([0 4])
    lgd = legend([hE_LQ, hE_MLQ_R, hE_MLQ_A, hE_MLQ_R2, hE_MLQ_A2, ...
        hE_MLQ_DR, hE_MLQ_DA], ...
        'LQ (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaA)})', ...
        'Location', 'best');
    title(lgd, [pattern{j} ' irradiaton'])
    grid on
    set(gca, 'FontSize', 16)

    h_MAE = figure();
    set(h_MAE, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    hE_LQ     = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.MAE.avg, ...
        '-og', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'g');
    hold on
    hE_MLQ_R  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R.MAE.avg, ...
        '-ob', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'b');
    hE_MLQ_A  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A.MAE.avg, ...
        '-or', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'r');
    hE_MLQ_R2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R2.MAE.avg, ...
        '-oc', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'c');
    hE_MLQ_A2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A2.MAE.avg, ...
        '-om', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'm');
    hE_MLQ_DR = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DR.MAE.avg, ...
        '-oy', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'y');
    hE_MLQ_DA = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DA.MAE.avg, ...
        '-ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    hold off
    xlabel('Dose (Gy)');
    ylabel('Mean Absolute Error');
    xlim([min(unique(dose_assay))-1 max(unique(dose_assay))+1])
    %     ylim([0 4])
    lgd = legend([hE_LQ, hE_MLQ_R, hE_MLQ_A, hE_MLQ_R2, hE_MLQ_A2, ...
        hE_MLQ_DR, hE_MLQ_DA], ...
        'LQ (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaA)})', ...
        'Location', 'best');
    title(lgd, [pattern{j} ' irradiaton'])
    grid on
    set(gca, 'FontSize', 16)

    h_RMSE = figure();
    set(h_RMSE, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    hE_LQ     = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.RMSE.avg, ...
        '-og', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'g');
    hold on
    hE_MLQ_R  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R.RMSE.avg, ...
        '-ob', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'b');
    hE_MLQ_A  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A.RMSE.avg, ...
        '-or', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'r');
    hE_MLQ_R2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R2.RMSE.avg, ...
        '-oc', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'c');
    hE_MLQ_A2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A2.RMSE.avg, ...
        '-om', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'm');
    hE_MLQ_DR = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DR.RMSE.avg, ...
        '-oy', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'y');
    hE_MLQ_DA = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DA.RMSE.avg, ...
        '-ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    hold off
    xlabel('Dose (Gy)');
    ylabel('Root-Mean-Squared Error');
    xlim([min(unique(dose_assay))-1 max(unique(dose_assay))+1])
    %     ylim([0 4])
    lgd = legend([hE_LQ, hE_MLQ_R, hE_MLQ_A, hE_MLQ_R2, hE_MLQ_A2, ...
        hE_MLQ_DR, hE_MLQ_DA], ...
        'LQ (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaA)})', ...
        'Location', 'best');
    title(lgd, [pattern{j} ' irradiaton'])
    grid on
    set(gca, 'FontSize', 16)

    h_LogLik = figure();
    set(h_LogLik, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    hE_LQ     = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.LogLik.avg, ...
        '-og', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'g');
    hold on
    hE_MLQ_R  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R.LogLik.avg, ...
        '-ob', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'b');
    hE_MLQ_A  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A.LogLik.avg, ...
        '-or', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'r');
    hE_MLQ_R2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R2.LogLik.avg, ...
        '-oc', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'c');
    hE_MLQ_A2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A2.LogLik.avg, ...
        '-om', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'm');
    hE_MLQ_DR = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DR.LogLik.avg, ...
        '-oy', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'y');
    hE_MLQ_DA = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DA.LogLik.avg, ...
        '-ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    hold off
    xlabel('Dose (Gy)');
    ylabel('Log-Likelihood');
    xlim([min(unique(dose_assay))-1 max(unique(dose_assay))+1])
    %     ylim([0 4])
    lgd = legend([hE_LQ, hE_MLQ_R, hE_MLQ_A, hE_MLQ_R2, hE_MLQ_A2, ...
        hE_MLQ_DR, hE_MLQ_DA], ...
        'LQ (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaA)})', ...
        'Location', 'best');
    title(lgd, [pattern{j} ' irradiaton'])
    grid on
    set(gca, 'FontSize', 16)

    h_Rsq = figure();
    set(h_Rsq, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    hE_LQ     = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).LQ.Rsq.avg, ...
        '-og', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'g');
    hold on
    hE_MLQ_R  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R.Rsq.avg, ...
        '-ob', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'b');
    hE_MLQ_A  = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A.Rsq.avg, ...
        '-or', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'r');
    hE_MLQ_R2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_R2.Rsq.avg, ...
        '-oc', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'c');
    hE_MLQ_A2 = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_A2.Rsq.avg, ...
        '-om', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'm');
    hE_MLQ_DR = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DR.Rsq.avg, ...
        '-oy', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'y');
    hE_MLQ_DA = plot(unique(dose_assay), PV_struct.(assay_pattern{j}).MLQ_DA.Rsq.avg, ...
        '-ok', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    hold off
    xlabel('Dose (Gy)');
    ylabel('R-squared');
    xlim([min(unique(dose_assay))-1 max(unique(dose_assay))+1])
%     ylim([0 1])
    lgd = legend([hE_LQ, hE_MLQ_R, hE_MLQ_A, hE_MLQ_R2, hE_MLQ_A2, ...
        hE_MLQ_DR, hE_MLQ_DA], ...
        'LQ (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaA^2)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaR)})', ...
        'MLQ (\it{\lambda_0exp(-\alphaD+\etaA)})', ...
        'Location', 'best');
    title(lgd, [pattern{j} ' irradiaton'])
    grid on
    set(gca, 'FontSize', 16)

    saveas(h_pred, fullfile(path, sprintf('%s_pred.png', assay_pattern{j})))
    saveas(h_MAE, fullfile(path, sprintf('%s_MAE.png', assay_pattern{j})))
    saveas(h_RMSE, fullfile(path, sprintf('%s_RMSE.png', assay_pattern{j})))
    saveas(h_LogLik, fullfile(path, sprintf('%s_LogLik.png', assay_pattern{j})))
    saveas(h_Rsq, fullfile(path, sprintf('%s_Rsq.png', assay_pattern{j})))

end

% plot(dose_vec.Control, count_vec.Control, 'ksquare', ...
%     dose_vec.Open, count_vec.Open, 'ro', ...
%     dose_vec.GRID_Stripes, count_vec.GRID_Stripes, 'g+', ...
%     dose_vec.GRID_Dots, count_vec.GRID_Dots, 'b^', ...
%     X_D, mdl_LQ.Fitted{:,1}, 'k-', ...
%     X_D, mdl_LQ_R.Fitted{:,1}, 'r-')
% hold on
%
% hold off
% xlabel('Dose (Gy)')
% ylabel('Colony count per quadrat; \lambda')
% lgd = legend({'Control', 'Open', 'Striped' 'Dotted', ...
%     'LQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
%     'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})'});
% title(lgd, 'Colony survival data')
% grid on
% set(gca, 'FontSize', 16)

%% Linear-quadratic (LQ) regression

replicates  = length(dose_assay)/length(unique(dose_assay));
D           = zeros(1, replicates);
SF          = ones(1, replicates);

for K = 1:length(dose_assay)

    D_temp      = dose_assay(ind_sort(K));
    SF_temp     = (assay.Open.img{ind_sort(K)}.N_colonies * f(D_temp)) ...
        / (N_colonies_ctrl * f(0)) ;  % (cells_seeded * PE_ctrl)]; %

    %     assay.open.img{ind_sort(k)}.SF = SF_temp;
    SF      = [SF SF_temp];
    D       = [D D_temp];

end

% LQ regression fit
tbl_LQ  = table(D', log(SF)', 'VariableNames', {'D', 'SF'});
mdl_LQ  = fitlm(tbl_LQ, 'SF ~ D + D^2', 'Intercept', false);

% Store alpha and beta LQ-parameters
alpha   = abs(mdl_LQ.Coefficients.Estimate(1));
beta    = abs(mdl_LQ.Coefficients.Estimate(2));
D_fit   = linspace(0, 10, 1001);
SF_fit  = exp(-alpha.*D_fit - beta.*D_fit.^2);

figure();
semilogy(D, SF, 'o', D_fit, SF_fit, '-r')
xlabel('Dose (Gy)')
ylabel('Surviving fraction')
legend('Dose response data', 'LQ fit (\it{exp(-\alphaD-\beta D^2)})', ...
    'Location', 'southwest')
grid on
set(gca, 'FontSize', 16)

%%  Plot EBT3 2D surviving fraction (SF) maps

% for K = 1:length(unique(dose_assay))
%     for i = 1:length(pattern)
%
%         % Estimate predicted survival
%         struct.(pattern{i}).(assay_dose{K}).img_SF = LQ( ...
%             dose_scale(K).*struct.(pattern{i}).img_dose, alpha, beta);
%
%         % Plot predicted survival map
%         plot2Dmap(struct.(pattern{i}).(assay_dose{K}).img_SF, ...
%             px_size, [0 1], 'Surviving fraction')
%
%     end
% end

%% Estimate predicted LQ survival profiles

% for K = 1:length(unique(dose_assay))
%
%     h_profilevec    = [];
%
%     for i = 1:length(pattern)
%
%         SFprofiles = [];
%         if strcmp('Open', pattern{i})
%             SF_EBT3.(pattern{i}).SFvalues        = [];
%         else
%             SF_EBT3.(pattern{i}).valley.SFvalues = [];
%             SF_EBT3.(pattern{i}).peak.SFvalues   = [];
%         end
%
%         for j = 1:size(struct.(pattern{i}).(assay_dose{K}).img_SF, 3)
%
%             % SF profiles
%             temp_img = struct.(pattern{i}).(assay_dose{K}). ...
%                 img_SF(1:end, x1_EBT3:x2_EBT3, j);
%             temp_img(temp_img <= 0) = 0;
%             temp_img(temp_img >= 1) = 1;
%
%             temp_SFprofile  = mean(temp_img, 2);
%             position        = (1:length(temp_SFprofile)) .* px_size;
%
%             % Save peak and valley SF for each profiles of a pattern
%             if strcmp('Open', pattern{i})
%                 SF_EBT3.(pattern{i}).SFvalues = [ ...
%                     SF_EBT3.(pattern{i}).SFvalues ...
%                     mean(temp_SFprofile(pos.(pattern{i})))];
%             else
%                 SF_EBT3.(pattern{i}).valley.SFvalues = [ ...
%                     SF_EBT3.(pattern{i}).valley.SFvalues ...
%                     mean(temp_SFprofile(pos.(pattern{i}).valley))];
%                 SF_EBT3.(pattern{i}).peak.SFvalues = [ ...
%                     SF_EBT3.(pattern{i}).peak.SFvalues ...
%                     mean(temp_SFprofile(pos.(pattern{i}).peak))];
%             end
%
%             % Store SF profile
%             SFprofiles = [SFprofiles; temp_SFprofile'];
%
%         end
%
%         if strcmp('Open', pattern{i})
%             [SF_EBT3.(pattern{i}).SFavg, SF_EBT3.(pattern{i}).CI] = ...
%                 estimateCI(SF_EBT3.(pattern{i}).SFvalues, significance);
%         else
%             [SF_EBT3.(pattern{i}).valley.SFavg, ...
%                 SF_EBT3.(pattern{i}).valley.CI] = ...
%                 estimateCI(SF_EBT3.(pattern{i}).valley.SFvalues, significance);
%             [SF_EBT3.(pattern{i}).peak.SFavg, ...
%                 SF_EBT3.(pattern{i}).peak.CI] = ...
%                 estimateCI(SF_EBT3.(pattern{i}).peak.SFvalues, significance);
%         end
%
%         % Estimate 95% CI bands for the SF profiles
%         profile.(pattern{i}).SF_LQ.(assay_dose{K}) = ...
%             estimateCIProfileBand(SFprofiles, position, significance);
%
%     end
%
% end

%% Estimate 1D survival profiles in bands

% clc
% shift = [0 -0.13 -0.15 0]; % [0 0.5 -0.075 0]; % -0.15 Stripes,
% temppos.GRID_Stripes.valley = [0.45 1.46 2.02 3.05 3.59 4.63 5.19];
% temppos.GRID_Dots.valley    = [0.51 1.08 3.15 3.70 5.80];
%
% for K = 1:length(unique(dose_assay))
%
%     for i = 2:length(assay_pattern)
%
%         counter = i - 2;
%         profile_count_vec   = [];
%         temp_posvec = [];
%         h_vec = [];
%
%         for j = 1:replicates
%
%             %             assay.(assay_pattern{i}).img{ind_sort((k-1)*replicates+j)}.filename
%             %             assay.(assay_pattern{i}).img{ind_sort((k-1)*replicates+j)}.dose
%             position = shift(i) + ...
%                 assay.(assay_pattern{i}).img{ind_sort((K-1)*replicates+j)}.position_band;
%             profile_count_vec = [profile_count_vec; ...
%                 assay.(assay_pattern{i}).img{ind_sort((K-1)*replicates+j)}.colony_count_band];
%
%         end
%
%         % Estimate 95% CI bands for the SF profiles, per dose
%         profile.(assay_pattern{i}).SF_assay.(assay_dose{K}) = ...
%             estimateCIProfileBand(profile_count_vec, position, 0.5);
%
%         % Normalization; estimate SF
%         if strcmp('Open', assay_pattern{i})
%             [~, ind_temp] = ...
%                 min(abs(D_calib-dose_scale(K)*round(dose_EBT3.(pattern{end-counter}).doseavg,2)));
%             profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).avgprofile = ...
%                 (profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).avgprofile ...
%                 .* f(ind_temp)) ./ (avg_colonycount_ctrl_band .* f(1));
%             profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_value = ...
%                 (profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_value ...
%                 .* f(ind_temp)) ./ (avg_colonycount_ctrl_band .* f(1));
%         else
%             temp_pos = ...
%                 profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).position;
%             temp_pos_CI = ...
%                 profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_pos;
%             if strcmp('GRID_Stripes', assay_pattern{i})
%                 ind_valley = ...
%                     (temp_pos > temppos.(assay_pattern{i}).valley(1) & ...
%                     temp_pos < temppos.(assay_pattern{i}).valley(2)) | ...
%                     (temp_pos > temppos.(assay_pattern{i}).valley(3) & ...
%                     temp_pos < temppos.(assay_pattern{i}).valley(4)) | ...
%                     (temp_pos > temppos.(assay_pattern{i}).valley(5) & ...
%                     temp_pos < temppos.(assay_pattern{i}).valley(6)) | ...
%                     (temp_pos > temppos.(assay_pattern{i}).valley(7));
%                 ind_valley_CI = ...
%                     (temp_pos_CI > temppos.(assay_pattern{i}).valley(1) & ...
%                     temp_pos_CI < temppos.(assay_pattern{i}).valley(2)) | ...
%                     (temp_pos_CI > temppos.(assay_pattern{i}).valley(3) & ...
%                     temp_pos_CI < temppos.(assay_pattern{i}).valley(4)) | ...
%                     (temp_pos_CI > temppos.(assay_pattern{i}).valley(5) & ...
%                     temp_pos_CI < temppos.(assay_pattern{i}).valley(6)) | ...
%                     (temp_pos_CI > temppos.(assay_pattern{i}).valley(7));
%             elseif strcmp('GRID_Dots', assay_pattern{i})
%                 ind_valley = ...
%                     (temp_pos < temppos.(assay_pattern{i}).valley(1) | ...
%                     temp_pos > temppos.(assay_pattern{i}).valley(2)) & ...
%                     (temp_pos < temppos.(assay_pattern{i}).valley(3) | ...
%                     temp_pos > temppos.(assay_pattern{i}).valley(4)) & ...
%                     (temp_pos < temppos.(assay_pattern{i}).valley(5));
%                 ind_valley_CI = ...
%                     (temp_pos_CI < temppos.(assay_pattern{i}).valley(1) | ...
%                     temp_pos_CI > temppos.(assay_pattern{i}).valley(2)) & ...
%                     (temp_pos_CI < temppos.(assay_pattern{i}).valley(3) | ...
%                     temp_pos_CI > temppos.(assay_pattern{i}).valley(4)) & ...
%                     (temp_pos_CI < temppos.(assay_pattern{i}).valley(5));
%             end
%             ind_peak = ~ind_valley;
%             ind_peak_CI = ~ind_valley_CI;
%
%             [~, ind_temp_valleydose] = ...
%                 min(abs(D_calib-dose_scale(K)*round(dose_EBT3.(pattern{end-counter}).valley.doseavg,2)));
%             [~, ind_temp_peakdose]   = ...
%                 min(abs(D_calib-dose_scale(K)*round(dose_EBT3.(pattern{end-counter}).peak.doseavg,2)));
%
%             % Valley
%             profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).avgprofile(ind_valley) = ...
%                 (profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).avgprofile(ind_valley) ...
%                 .* f(ind_temp_valleydose)) ./ (avg_colonycount_ctrl_band .* f(1));
%             profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_value(ind_valley_CI) = ...
%                 (profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_value(ind_valley_CI) ...
%                 .* f(ind_temp_valleydose)) ./ (avg_colonycount_ctrl_band .* f(1));
%
%             % Peak
%             profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).avgprofile(ind_peak) = ...
%                 (profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).avgprofile(ind_peak) ...
%                 .* f(ind_temp_peakdose)) ./ (avg_colonycount_ctrl_band .* f(1));
%             profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_value(ind_peak_CI) = ...
%                 (profile.(assay_pattern{i}).SF_assay.(assay_dose{K}).CI95_value(ind_peak_CI) ...
%                 .* f(ind_temp_peakdose)) ./ (avg_colonycount_ctrl_band .* f(1));
%
%         end
%
%     end
%
% end

% for k = 1:length(unique(dose_assay))
%     for i = 3 % 1:length(pattern)
%
%         h_vec = [];
%
%         h = figure();
%         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%         hold on
%         h_95CIband_observed = fill( ...
%             profile.(pattern{i}).SF_assay.(assay_dose{k}).CI95_pos, ...
%             profile.(pattern{i}).SF_assay.(assay_dose{k}).CI95_value, ...
%             [0.5 0.5 0.5], 'FaceColor', [0.7 0.7 0.7], ...
%             'EdgeColor', 'none', 'FaceAlpha', 0.5);
%         h_95CImean_observed = plot( ...
%             profile.(pattern{i}).SF_assay.(assay_dose{k}).position, ...
%             profile.(pattern{i}).SF_assay.(assay_dose{k}).avgprofile, ...
%             '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
%         ylabel('Surviving fraction')
%         xlabel('Position (cm)')
%         xlim([0.25 5.75])
%         ylim(y_limits_SF{k})
%         h_vec = [h_vec h_95CImean_observed h_95CIband_observed];
%         lgd = legend(h_vec, 'Observed survival', 'Location', 'NorthEast', ...
%             'NumColumns', 3, 'Interpreter', 'LaTeX');
%         title(lgd, [strrep(pattern{i}, '_', ' ') ' irradiaton, ', ...
%             sprintf('%.1f',  dose_scale(k)*5) ' Gy nominal dose'] )
%                 grid on
%         set(gca, 'FontSize', 26)
%         hold off
%
%     end
% end

%% Plot 95% CI of SF profiles for observed vs. predicted

% for K = 1:length(unique(dose_assay))
%
%     for i = 1:length(pattern)
%
%         h_vec = [];
%
%         temp_SF_LQ = interp1( ...
%             profile.(pattern{i}).SF_LQ.(assay_dose{K}).position, ...
%             profile.(pattern{i}).SF_LQ.(assay_dose{K}).avgprofile, ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).position);
%         RPD = 100 .* (temp_SF_LQ - ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).avgprofile) ./ ...
%             ( (profile.(pattern{i}).SF_assay.(assay_dose{K}).avgprofile + ...
%             temp_SF_LQ)./2 );
%
%         if strcmp('Open', pattern{i})
%             temp_pos = ...
%                 profile.(pattern{i}).SF_assay.(assay_dose{K}).position;
%             ind_Open = temp_pos > 0.5 & temp_pos < 5.1;
%             pattern{i}
%             mean(RPD(ind_Open))
%         else
%             temp_pos = ...
%                 profile.(pattern{i}).SF_assay.(assay_dose{K}).position;
%             ind_valley = (temp_pos > 0.55 & temp_pos < 1.35) | ...
%                 (temp_pos > 2.1 & temp_pos < 2.95) | ...
%                 (temp_pos > 3.7 & temp_pos < 4.55);
%             ind_peak = (temp_pos > 1.55 & temp_pos < 1.95) | ...
%                 (temp_pos > 3.12 & temp_pos < 3.5) | ...
%                 (temp_pos > 4.7 & temp_pos < 5.1);
%             pattern{i}
%             mean(RPD(ind_valley))
%             mean(RPD(ind_peak))
%         end
%
%         h = figure();
%         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%         hold on
%
%         yyaxis left
%         % Plot predicted response
%         h_95CIband_predicted = fill( ...
%             profile.(pattern{i}).SF_LQ.(assay_dose{K}).CI95_pos, ...
%             profile.(pattern{i}).SF_LQ.(assay_dose{K}).CI95_value, ...
%             color_CI95{i}, 'FaceColor', facecolor_CI95{i}, ...
%             'EdgeColor', 'none');
%         h_95CImean_predicted = plot( ...
%             profile.(pattern{i}).SF_LQ.(assay_dose{K}).position, ...
%             profile.(pattern{i}).SF_LQ.(assay_dose{K}).avgprofile, ...
%             linespec{i}, 'LineWidth', 1.0);
%         % Plot observed response
%         h_95CIband_observed = fill( ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).CI95_pos, ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).CI95_value, ...
%             [0.5 0.5 0.5], 'FaceColor', [0.7 0.7 0.7], ...
%             'EdgeColor', 'none', 'FaceAlpha', 0.5);
%         h_95CImean_observed = plot( ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).position, ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).avgprofile, ...
%             '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
%         ylabel('Surviving fraction')
%         ylim(y_limits_SF{K})
%         set(gca, 'FontSize', 26) %, 'YScale', 'log')
%
%         % Plot relative percentage difference
%         yyaxis right
%         h_RPE = plot( ...
%             profile.(pattern{i}).SF_assay.(assay_dose{K}).position, ...
%             RPD, '--k');
%         ylabel('Relative percentage deviation (%)')
%         ylim(y_limits_RPE{K})
%         yline(0, '-', 'RPD = 0', 'LineWidth', 1, 'FontSize', 18);
%         h_vec = [h_vec h_95CImean_predicted h_95CIband_predicted ...
%             h_95CImean_observed h_95CIband_observed h_RPE];
%
%         ax = gca;
%         ax.YAxis(1).Color = 'k';
%         ax.YAxis(2).Color = 'k';
%
%         xlabel('Position (cm)')
%         xlim([0.25 5.75])
%         lgd = legend(h_vec, lgdstr, 'Location', 'NorthEast', ...
%             'NumColumns', 3, 'Interpreter', 'LaTeX');
%         title(lgd, [strrep(pattern{i}, '_', ' ') ' irradiaton, ', ...
%             sprintf('%.1f',  dose_scale(K)*5) ' Gy nominal dose'] )
%         %         grid on
%         set(gca, 'FontSize', 26)
%         hold off
%
%     end
%
% end
