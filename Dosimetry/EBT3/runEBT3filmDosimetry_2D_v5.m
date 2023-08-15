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
path_assay_img      = 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549';
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
irradarea.Control      = 0.0;
irradarea.GRID_Dots    = 0.05;
irradarea.GRID_Stripes = 0.33;
irradarea.Open         = 1.0;

peakdose.Control      = 0.0;        valleydose.Control      = 0.0;
peakdose.GRID_Dots    = 0.70;       valleydose.GRID_Dots    = 0.10;
peakdose.GRID_Stripes = 0.82;       valleydose.GRID_Stripes = 0.18;
peakdose.Open         = 1.0;        valleydose.Open         = 0.0;

for j = 1:length(pattern)

    % Store dose map
    temp_img_dose = mean(struct.(pattern{j}).img_dose, 3);
    temp_img_dose(temp_img_dose <= 0)       = 0;
    temp_img_dose(isnan(temp_img_dose))     = 0;
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
    if strcmp(pattern{j}, 'GRID_Dots')
        
        % Valley-to-Peak distance map
        bw_temp = zeros(size(temp_bw_dose,1), size(temp_bw_dose,2));
        for i = 1:size(peak_xycoord.(pattern{j}),1)
            bw_temp(peak_xycoord.(pattern{j})(i,2), ...
                    peak_xycoord.(pattern{j})(i,1)) = 1;
        end
        img_valleytopeak_dist = bwdist(bw_temp);
%         figure(); imshow(bw_temp, [])

        % Peak-to-Valley distance map
        bw_temp = zeros(size(temp_bw_dose,1), size(temp_bw_dose,2));
        [colsInImage, rowsInImage] = meshgrid(1:size(temp_bw_dose,2), ...
                                              1:size(temp_bw_dose,1));
        for i = 1:size(peak_xycoord.(pattern{j}),1)
            circlePixels = ...
                (rowsInImage - peak_xycoord.(pattern{j})(i,2)).^2 + ...
                (colsInImage - peak_xycoord.(pattern{j})(i,1)).^2 <= 105^2;
            bw_temp = logical(bw_temp | circlePixels);
        end
        img_peaktovalley_dist = bwdist(~bw_temp);
%         figure(); imshow(bw_temp, [])

    elseif strcmp(pattern{j}, 'GRID_Stripes')

        % Valley-to-Peak distance map
        bw_temp = zeros(size(temp_bw_dose,1), size(temp_bw_dose,2));
        bw_temp(peak_xycoord.(pattern{j})(:,2), ...
                peak_xycoord.(pattern{j})(:,1)) = 1;
        img_valleytopeak_dist = bwdist(bw_temp);
% %         figure(); imshow(bw_temp, [])

        % Peak-to-Valley distance map
        bw_temp = zeros(size(temp_bw_dose,1), size(temp_bw_dose,2));
        [colsInImage, rowsInImage] = meshgrid(1:size(temp_bw_dose,2), ...
                                              1:size(temp_bw_dose,1));
        for i = 1:size(peak_xycoord.(pattern{j}),1)
            circlePixels = ...
                (rowsInImage - peak_xycoord.(pattern{j})(i,2)).^2 + ...
                (colsInImage - peak_xycoord.(pattern{j})(i,1)).^2 <= 80^2;
            bw_temp = logical(bw_temp | circlePixels);
        end
        img_peaktovalley_dist = bwdist(~bw_temp);
%         figure(); imshow(bw_temp, [])

    else

        img_valleytopeak_dist = zeros(size(temp_img_dose,1), ...
                                      size(temp_img_dose,2));
        img_peaktovalley_dist = zeros(size(temp_img_dose,1), ...
                                      size(temp_img_dose,2));
        
    end
    img_valleytopeak_dist = px_size .* img_valleytopeak_dist .* temp_bw_dose;
    img_peaktovalley_dist = px_size .* img_peaktovalley_dist .* temp_bw_dose;
    
    % Compute peak area image
    img_irradarea = irradarea.(pattern{j}) .* ones(size(temp_img_dose,1), ...
        size(temp_img_dose,2));
    img_irradarea = temp_bw_dose .* img_irradarea;

    % Compute peak and valley dose image
    img_peakdose = peakdose.(pattern{j}) .* ones(size(temp_img_dose,1), ...
        size(temp_img_dose,2));
    img_valleydose = valleydose.(pattern{j}) .* ones(size(temp_img_dose,1), ...
        size(temp_img_dose,2));
    img_peakdose    = temp_bw_dose .* img_peakdose;
    img_valleydose  = temp_bw_dose .* img_valleydose;

    % Store distance, area and dose images
    struct.(pattern{j}).img_valleytopeak_dist = img_valleytopeak_dist;
    struct.(pattern{j}).img_peaktovalley_dist = img_peaktovalley_dist;
    struct.(pattern{j}).img_irradarea         = img_irradarea;
    struct.(pattern{j}).img_peakdose          = img_peakdose;
    struct.(pattern{j}).img_valleydose        = img_valleydose;
    %     struct.(pattern{j}).img_grad     = img_grad;

    % Plot EBT3 2D dose and peak distance map
    plot2Dmap(struct.(pattern{j}).img_dose, px_size, [0 6], ...
        'Dose (Gy)')
    plot2Dmap(struct.(pattern{j}).img_peaktovalley_dist, px_size, ...
        [0 max(max(struct.(pattern{j}).img_peaktovalley_dist))], ...
        'Peak-to-Valley distance (cm)')
    plot2Dmap(struct.(pattern{j}).img_valleytopeak_dist, px_size, ...
        [0 max(max(struct.(pattern{j}).img_valleytopeak_dist))], ...
        'Valley-to-Peak distance (cm)')
    plot2Dmap(img_irradarea, px_size, ...
        [0 max(max(struct.(pattern{j}).img_irradarea))], ...
        'Irradiation area fraction')
    %     plot2Dmap(struct.(pattern{j}).img_grad, px_size, ...
    %         [0 max(max(struct.(pattern{j}).img_grad))], ...
    %         'Dose gradient (cm/Gy)')

end

%%

path = 'C:\Users\delmo\Desktop\New folder1\';

px_size_assay   = 2.54/1200;    % in cm/px, 1200 dpi
cells_seeded    = 10000;        % no. of seeded cells (30000 originally)

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
translation_vec_seg.Control   = {[2 19], [4 31], [0 0], [2 39]};
translation_vec_seg.GRID_Dots = {[-12 -4], [0 0], [-5 22], [1 5], ...
    [-9 -1], [-11 22], [-9 -4], [-3 2], ...
    [-5 31], [-3 31], [-5 31], [6 34]};
translation_vec_seg.GRID_Stripes = {[-5 28], [-2 28], [-5 15], [0 0], ...
    [-4 25], [-5 28], [-3 44], [-3 35], ...
    [-2 18], [-4 -1], [-2 2], [-6 30]};
translation_vec_seg.Open = {[-4 35], [3 25], [0 20], [-2 2], ...
    [5 51], [2 42], [-4 3], [-6 14], ...
    [-9 43], [-3 44], [-6 -7], [2 3]};

translation_vec.GRID_Dots = {[-116 55], [5 78], [-34 71], [-38 -1], ...
    [-70 45], [10 48], [85 8], [-16 1], ...
    [-50 43], [0 59], [89 -5], [20 -15]};
translation_vec.GRID_Stripes = {[0 67], [0 60], [0 -15], [0 55], ...
    [0 25], [0 35], [0 -15], [0 -25], ...
    [0 41], [0 48], [0 11], [0 -20]};

translation_dose_vec.GRID_Dots = [-41 37];
translation_dose_vec.GRID_Stripes = [96 85];

folderlist_assay     = getAllFolders(path_assay);
folderlist_assay_img = getAllFolders(path_assay_img);
folderlist_bwflask   = getAllFolders(path_bwflask);

close all
clc
for i = 1:length(folderlist_assay)

    [~, pattern_name, ~] = fileparts(folderlist_assay{i});
    filelist_assay       = getAllFiles(folderlist_assay{i});
    filelist_assay_img   = getAllFiles(folderlist_assay_img{i});
    filelist_bwflask     = getAllFiles(folderlist_bwflask{i});

    if i > 1
        temp_img_dose = mean(struct.(assay_pattern{i}).img_dose, 3);
        temp_img_dose(temp_img_dose <= 0) = 0;
        temp_img_dose(isnan(temp_img_dose)) = 0;
        temp_img_dose(imag(temp_img_dose) ~= 0) = 0;
        temp_img_dose = medfilt2(temp_img_dose, [5 5]);
        temp_img_peaktovalley_dist = struct.(assay_pattern{i}).img_peaktovalley_dist;
        temp_img_valleytopeak_dist = struct.(assay_pattern{i}).img_valleytopeak_dist;
        temp_img_irradarea         = struct.(assay_pattern{i}).img_irradarea;
        temp_img_peakdose          = struct.(assay_pattern{i}).img_peakdose;
        temp_img_valleydose        = struct.(assay_pattern{i}).img_valleydose;
    end

    for j = 1:length(filelist_assay)

        [~, fn_assay, ~]     = fileparts(filelist_assay{j});
        [~, fn_assay_img, ~] = fileparts(filelist_assay_img{j});
        bw_assay_seg         = logical(readmatrix(filelist_assay{j}));
        bw_flask             = logical(readmatrix(filelist_bwflask{j}));
        [img_assay, ~]       = imread(filelist_assay_img{j});
        img_assay            = img_assay(:,:,1:3);

        % Manual image registration for
        bw_assay_seg = imtranslate(bw_assay_seg, ...
            translation_vec_seg.(assay_pattern{i}){j});

%         h = figure();
%         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%         imshow(img_assay)
%         hold on
%         for row = 1:dxdy_px(1):3008 % size(assay.(assay_pattern{i}).img{j}.img_dose,1)
%             line([1, 2209], [row, row], 'Color', [.7 .7 .7])
%         end
%         for col = 1:dxdy_px(2):2209 % size(assay.(assay_pattern{i}).img{j}.img_dose,2)
%             line([col, col], [1, 3008], 'Color', [.7 .7 .7])
%         end
%         visboundaries(bw_flask, 'Color', 'yellow', 'LineWidth', 1)
%         visboundaries(bw_assay_seg, 'Color', 'red', 'LineWidth', 1)
%         hold off
%         set(gca, 'FontSize', 16)
        %         saveas(h, fullfile(path, [fn_assay_img '_imgassay.tiff']))

        % Manual image registration
        if contains(pattern_name, ["Dots", "Stripes"]) % && j > 4
            bw_assay_seg = imtranslate(bw_assay_seg, translation_vec.(assay_pattern{i}){j});
            bw_flask     = imtranslate(bw_flask, translation_vec.(assay_pattern{i}){j});
            img_assay    = imtranslate(img_assay, translation_vec.(assay_pattern{i}){j});
        end

        stats       = regionprops(bw_assay_seg, 'Centroid');
        centroids   = round(cat(1, stats.Centroid));
        x_centroids = centroids(:,2);
        y_centroids = centroids(:,1);
        N_colonies  = length(stats);
        assay.(assay_pattern{i}).img{j}.filename    = fn_assay_img;
        assay.(assay_pattern{i}).img{j}.bw_seg      = bw_assay_seg;
        assay.(assay_pattern{i}).img{j}.bw_flask    = bw_flask;
        assay.(assay_pattern{i}).img{j}.x_centroids = x_centroids;
        assay.(assay_pattern{i}).img{j}.y_centroids = y_centroids;
        assay.(assay_pattern{i}).img{j}.N_colonies_uncalib = length(stats);
        if strcmp(assay_pattern{i}, 'Control')
            assay.(assay_pattern{i}).img{j}.PE      = ...
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
        assay.(assay_pattern{i}).img{j}.bw_flask         = bw_flask;
        assay.(assay_pattern{i}).img{j}.bw_seg           = bw_assay_seg;
        assay.(assay_pattern{i}).img{j}.bw_centroids     = bw_centroids;
        assay.(assay_pattern{i}).img{j}.n_quadrat_x      = n_quadrat_x;
        assay.(assay_pattern{i}).img{j}.n_quadrat_y      = n_quadrat_y;
        assay.(assay_pattern{i}).img{j}.bw_flask_quadrat = temp_bw_flask_quadrat;
        assay.(assay_pattern{i}).img{j}.bw_seg_quadrat   = bw_assay_seg_quadrat;

        % Scale dose map accordingly to nominal dose
        if i == 1
            temp_img_dose_scaled = 0.0 .* temp_img_dose;
        else
            temp_img_dose_scaled       = (dose_assay(j)/5) .* temp_img_dose;
            temp_img_peakdose_scaled   = dose_assay(j) .* temp_img_peakdose;
            temp_img_valleydose_scaled = dose_assay(j) .* temp_img_valleydose;
        end

        % Interpolate dose maps to match spatial resolution of assay images
        [X,Y]   = meshgrid(1:size(temp_img_dose_scaled,2), 1:size(temp_img_dose_scaled,1));
        [X2,Y2] = meshgrid( ...
            1:px_size_assay/px_size:size(temp_img_dose_scaled,2), ...
            1:px_size_assay/px_size:size(temp_img_dose_scaled,1));
        interp_img_dose = interp2(X, Y, temp_img_dose_scaled, X2, Y2, 'spline');
        dummy_zeros_x = zeros(size(bw_assay_seg,1)-size(interp_img_dose,1), ...
            size(interp_img_dose,2));
        dummy_zeros_y = zeros(size(bw_assay_seg,1), ...
            size(bw_assay_seg,2)-size(interp_img_dose,2));
        interp_img_dose = cat(1, interp_img_dose, dummy_zeros_x);
        interp_img_dose = cat(2, interp_img_dose, dummy_zeros_y);
        interp_img_dose(interp_img_dose < 0) = 0;

        if i == 1

            img_quadrat_zeros      = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_count_quadrat = round(f(0) .* temp_img_count_quadrat);
            temp_img_count_quadrat = temp_img_count_quadrat .* temp_bw_flask_quadrat;

            % Store quadrat results
            assay.(assay_pattern{i}).img{j}.dose_avg               = 0;
            assay.(assay_pattern{i}).img{j}.dose_nom               = 0;
            assay.(assay_pattern{i}).img{j}.N_colonies             = sum(temp_img_count_quadrat(:)); % f(0)*N_colonies;
            assay.(assay_pattern{i}).img{j}.N_quadrat              = nnz(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.bw_quadrat             = temp_bw_flask_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose_quadrat       = img_quadrat_zeros;
            assay.(assay_pattern{i}).img{j}.img_peaktovalley_dist_quadrat = img_quadrat_zeros;
            assay.(assay_pattern{i}).img{j}.img_valleytopeak_dist_quadrat = img_quadrat_zeros;
            assay.(assay_pattern{i}).img{j}.img_irradarea_quadrat  = img_quadrat_zeros;
            assay.(assay_pattern{i}).img{j}.img_peakdose_quadrat   = img_quadrat_zeros;
            assay.(assay_pattern{i}).img{j}.img_valleydose_quadrat = img_quadrat_zeros;
            assay.(assay_pattern{i}).img{j}.img_count_quadrat      = temp_img_count_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose               = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_peaktovalley_dist  = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_valleytopeak_dist  = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_irradarea          = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_peakdose           = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_valleydose         = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.PV.dose_quadrat        = img_quadrat_zeros(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.peaktovalley_dist_quadrat = img_quadrat_zeros(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.valleytopeak_dist_quadrat = img_quadrat_zeros(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.irradarea_quadrat   = img_quadrat_zeros(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.peakdose_quadrat    = img_quadrat_zeros(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.valleydose_quadrat  = img_quadrat_zeros(temp_bw_flask_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.count_quadrat       = temp_img_count_quadrat(temp_bw_flask_quadrat);

            % Plot quadrat images
            h = figure();
            set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
            tlo = tiledlayout(3,3);
            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_dose_quadrat, ...
                [0 dose_assay(j)], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_peaktovalley_dist_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Peak-to-Valley distance (cm)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_irradarea_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Irradiation fraction (a.u.)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_peakdose_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Peak dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_valleydose_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Valey dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_count_quadrat, ...
                [0 max(temp_img_count_quadrat(:))], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Surviving Colony Count (a.u.)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

        else

            % Interpolate distance and irradiation area maps
            interp_img_peaktovalley_dist = interp2(X, Y, temp_img_peaktovalley_dist, X2, Y2, 'spline');
            interp_img_peaktovalley_dist = cat(1, interp_img_peaktovalley_dist, dummy_zeros_x);
            interp_img_peaktovalley_dist = cat(2, interp_img_peaktovalley_dist, dummy_zeros_y);

            interp_img_valleytopeak_dist = interp2(X, Y, temp_img_valleytopeak_dist, X2, Y2, 'spline');
            interp_img_valleytopeak_dist = cat(1, interp_img_valleytopeak_dist, dummy_zeros_x);
            interp_img_valleytopeak_dist = cat(2, interp_img_valleytopeak_dist, dummy_zeros_y);

            interp_img_irradarea = interp2(X, Y, temp_img_irradarea, X2, Y2, 'spline');
            interp_img_irradarea = cat(1, interp_img_irradarea, dummy_zeros_x);
            interp_img_irradarea = cat(2, interp_img_irradarea, dummy_zeros_y);

            interp_img_peakdose = interp2(X, Y, temp_img_peakdose_scaled, X2, Y2, 'spline');
            interp_img_peakdose = cat(1, interp_img_peakdose, dummy_zeros_x);
            interp_img_peakdose = cat(2, interp_img_peakdose, dummy_zeros_y);

            interp_img_valleydose = interp2(X, Y, temp_img_valleydose_scaled, X2, Y2, 'spline');
            interp_img_valleydose = cat(1, interp_img_valleydose, dummy_zeros_x);
            interp_img_valleydose = cat(2, interp_img_valleydose, dummy_zeros_y);

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
                bw_dose = imwarp(bw_dose, tform, 'OutputView', Rfixed);
                interp_img_dose = imwarp(interp_img_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_peaktovalley_dist = imwarp(interp_img_peaktovalley_dist, tform, ...
                    'OutputView', Rfixed);
                interp_img_valleytopeak_dist = imwarp(interp_img_valleytopeak_dist, tform, ...
                    'OutputView', Rfixed);
                interp_img_irradarea = imwarp(interp_img_irradarea, tform, ...
                    'OutputView', Rfixed);
                interp_img_peakdose = imwarp(interp_img_peakdose, tform, ...
                    'OutputView', Rfixed);
                interp_img_valleydose = imwarp(interp_img_valleydose, tform, ...
                    'OutputView', Rfixed);

            elseif j > 1 && strcmp(assay_pattern{i}, 'Open')

                bw_dose = imwarp(bw_dose, tform, 'OutputView', Rfixed);
                interp_img_dose = imwarp(interp_img_dose, tform, ...
                    'OutputView', Rfixed);
                interp_img_peaktovalley_dist = imwarp(interp_img_peaktovalley_dist, tform, ...
                    'OutputView', Rfixed);
                interp_img_valleytopeak_dist = imwarp(interp_img_valleytopeak_dist, tform, ...
                    'OutputView', Rfixed);
                interp_img_irradarea = imwarp(interp_img_irradarea, tform, ...
                    'OutputView', Rfixed);
                interp_img_peakdose = imwarp(interp_img_peakdose, tform, ...
                    'OutputView', Rfixed);
                interp_img_valleydose = imwarp(interp_img_valleydose, tform, ...
                    'OutputView', Rfixed);

            elseif any(strcmp(assay_pattern{i}, {'GRID_Stripes', 'GRID_Dots'}))

                bw_dose = imtranslate(bw_dose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_dose = imtranslate(interp_img_dose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_peaktovalley_dist = imtranslate(interp_img_peaktovalley_dist, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_valleytopeak_dist = imtranslate(interp_img_valleytopeak_dist, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_irradarea = imtranslate(interp_img_irradarea, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_peakdose = imtranslate(interp_img_peakdose, ...
                    translation_dose_vec.(assay_pattern{i}));
                interp_img_valleydose = imtranslate(interp_img_valleydose, ...
                    translation_dose_vec.(assay_pattern{i}));

            end
            bw_dose = imerode(bw_dose, strel('disk', 20));
            bw_dose = bwmorph(bw_dose, 'shrink');
            bw_dose = bwmorph(bw_dose, 'thin');
            bw_dose = bwareafilt(bw_dose, 1);
            bw_dose = logical(bw_dose);
            %             figure(); imshow(bw_dose, [])

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
            img_peaktovalley_dist_quadrat = mat2quadrat(interp_img_peaktovalley_dist, dxdy_px);
            img_valleytopeak_dist_quadrat = mat2quadrat(interp_img_valleytopeak_dist, dxdy_px);
            img_irradarea_quadrat   = mat2quadrat(interp_img_irradarea, dxdy_px);
            img_peakdose_quadrat    = mat2quadrat(interp_img_peakdose, dxdy_px);
            img_valleydose_quadrat  = mat2quadrat(interp_img_valleydose, dxdy_px);

            temp_bw_dose_quadrat        = false(n_quadrat_x, n_quadrat_y);
            temp_bw_dose_peak_quadrat   = false(n_quadrat_x, n_quadrat_y);
            temp_bw_dose_valley_quadrat = false(n_quadrat_x, n_quadrat_y);
            temp_img_dose_quadrat       = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_peaktovalley_dist_quadrat = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_valleytopeak_dist_quadrat = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_irradarea_quadrat  = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_peakdose_quadrat   = zeros(n_quadrat_x, n_quadrat_y);
            temp_img_valleydose_quadrat = zeros(n_quadrat_x, n_quadrat_y);

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
                    temp_img_dose_quadrat(k,l)       = mean(img_dose_quadrat{k,l}(:));
                    temp_img_peaktovalley_dist_quadrat(k,l) = mean(img_peaktovalley_dist_quadrat{k,l}(:));
                    temp_img_valleytopeak_dist_quadrat(k,l) = mean(img_valleytopeak_dist_quadrat{k,l}(:));
                    temp_img_irradarea_quadrat(k,l)  = mean(img_irradarea_quadrat{k,l}(:));
                    temp_img_peakdose_quadrat(k,l)   = mean(img_peakdose_quadrat{k,l}(:));
                    temp_img_valleydose_quadrat(k,l) = mean(img_valleydose_quadrat{k,l}(:));
                end
            end
            temp_img_count_quadrat = round(f(temp_img_dose_quadrat) ...
                .* temp_img_count_quadrat);

            temp_bw                = bw_flask & bw_dose;
            temp_bw_peak           = temp_bw_dose_peak & temp_bw;
            temp_bw_valley         = temp_bw_dose_valley & temp_bw;
            temp_bw_quadrat        = temp_bw_flask_quadrat & temp_bw_dose_quadrat;
            temp_bw_peak_quadrat   = temp_bw_dose_peak_quadrat & temp_bw_quadrat;
            temp_bw_valley_quadrat = temp_bw_dose_valley_quadrat & temp_bw_quadrat;
            %            temp_bw_quadrat = temp_bw_dose_peak_quadrat | temp_bw_dose_valley_quadrat; % NEW!!!

            % Get relevant pixel and quadrat values
            temp_img_dose_quadrat       = temp_img_dose_quadrat .* temp_bw_quadrat;
            temp_img_peaktovalley_dist_quadrat = temp_img_peaktovalley_dist_quadrat .* temp_bw_quadrat;
            temp_img_valleytopeak_dist_quadrat = temp_img_valleytopeak_dist_quadrat .* temp_bw_quadrat;
            temp_img_irradarea_quadrat  = temp_img_irradarea_quadrat .* temp_bw_quadrat;
            temp_img_peakdose_quadrat   = temp_img_peakdose_quadrat .* temp_bw_quadrat;
            temp_img_valleydose_quadrat = temp_img_valleydose_quadrat .* temp_bw_quadrat;
            temp_img_count_quadrat      = temp_img_count_quadrat .* temp_bw_quadrat;

            % Store quadrat results
            assay.(assay_pattern{i}).img{j}.dose_avg               = mean(interp_img_dose(temp_bw));
            assay.(assay_pattern{i}).img{j}.dose_nom               = dose_assay(j);
            assay.(assay_pattern{i}).img{j}.N_colonies             = sum(temp_img_count_quadrat(:)); % f(mean(interp_img_dose(temp_bw))) * N_colonies;
            assay.(assay_pattern{i}).img{j}.N_quadrat              = nnz(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.bw_dose                = bw_dose;
            assay.(assay_pattern{i}).img{j}.bw_dose_quadrat        = temp_bw_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.bw_quadrat             = temp_bw_quadrat;
            assay.(assay_pattern{i}).img{j}.img_dose               = interp_img_dose;
            assay.(assay_pattern{i}).img{j}.img_peaktovalley_dist  = interp_img_peaktovalley_dist;
            assay.(assay_pattern{i}).img{j}.img_valleytopeak_dist  = interp_img_valleytopeak_dist;
            assay.(assay_pattern{i}).img{j}.img_irradarea          = interp_img_irradarea;
            assay.(assay_pattern{i}).img{j}.img_peakdose           = interp_img_peakdose;
            assay.(assay_pattern{i}).img{j}.img_valleydose         = interp_img_valleydose;
            assay.(assay_pattern{i}).img{j}.img_dose_quadrat       = temp_img_dose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_peaktovalley_dist_quadrat = temp_img_peaktovalley_dist_quadrat;
            assay.(assay_pattern{i}).img{j}.img_valleytopeak_dist_quadrat = temp_img_valleytopeak_dist_quadrat;
            assay.(assay_pattern{i}).img{j}.img_irradarea_quadrat  = temp_img_irradarea_quadrat;
            assay.(assay_pattern{i}).img{j}.img_peakdose_quadrat   = temp_img_peakdose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_valleydose_quadrat = temp_img_valleydose_quadrat;
            assay.(assay_pattern{i}).img{j}.img_count_quadrat      = temp_img_count_quadrat;
            assay.(assay_pattern{i}).img{j}.PV.dose_quadrat        = temp_img_dose_quadrat(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.peaktovalley_dist_quadrat = temp_img_peaktovalley_dist_quadrat(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.valleytopeak_dist_quadrat = temp_img_valleytopeak_dist_quadrat(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.irradarea_quadrat   = temp_img_irradarea_quadrat(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.peakdose_quadrat    = temp_img_peakdose_quadrat(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.valleydose_quadrat  = temp_img_valleydose_quadrat(temp_bw_quadrat);
            assay.(assay_pattern{i}).img{j}.PV.count_quadrat       = temp_img_count_quadrat(temp_bw_quadrat);
            if any(strcmp(assay_pattern{i}, {'GRID_Stripes', 'GRID_Dots'}))
                assay.(assay_pattern{i}).img{j}.N_quadrat_peak               = nnz(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.N_quadrat_valley             = nnz(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.bw_peak_quadrat              = temp_bw_peak_quadrat;
                assay.(assay_pattern{i}).img{j}.bw_valley_quadrat            = temp_bw_valley_quadrat;
                assay.(assay_pattern{i}).img{j}.PV.dose_peak_quadrat         = temp_img_dose_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.peaktovalley_dist_peak_quadrat = temp_img_peaktovalley_dist_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.valleytopeak_dist_peak_quadrat = temp_img_valleytopeak_dist_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.irradarea_peak_quadrat    = temp_img_irradarea_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.peakdose_peak_quadrat     = temp_img_peakdose_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.valleydose_peak_quadrat   = temp_img_valleydose_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.count_peak_quadrat        = temp_img_count_quadrat(temp_bw_peak_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.dose_valley_quadrat       = temp_img_dose_quadrat(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.peaktovalley_dist_valley_quadrat = temp_img_peaktovalley_dist_quadrat(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.valleytopeak_dist_valley_quadrat = temp_img_valleytopeak_dist_quadrat(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.irradarea_valley_quadrat  = temp_img_irradarea_quadrat(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.peakdose_valley_quadrat   = temp_img_peakdose_quadrat(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.valleydose_valley_quadrat = temp_img_valleydose_quadrat(temp_bw_valley_quadrat);
                assay.(assay_pattern{i}).img{j}.PV.count_valley_quadrat      = temp_img_count_quadrat(temp_bw_valley_quadrat);
            end

            % Plot quadrat images
            h = figure();
            set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
            tlo = tiledlayout(3,3);
            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_dose_quadrat, ...
                [0 dose_assay(j)], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

%             nexttile(tlo);
%             imshow(assay.(assay_pattern{i}).img{j}.img_peaktovalley_dist_quadrat, ...
%                 [], 'colormap', jet(4096))
%             shading interp
%             c = colorbar;
%             c.Label.String = 'Peak-to-Valley distance (cm)';
%             hold on
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
%                 'Color', 'yellow', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
%                 'Color', 'green', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
%                 'Color', 'red', 'LineWidth', 1)
%             hold off
%             set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_valleytopeak_dist_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Valley-to-Peak distance (cm)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_irradarea_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Irradiation fraction (a.u.)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_peakdose_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Peak dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_valleydose_quadrat, ...
                [], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Valey dose (Gy)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

            nexttile(tlo);
            imshow(assay.(assay_pattern{i}).img{j}.img_count_quadrat, ...
                [0 max(temp_img_count_quadrat(:))], 'colormap', jet(4096))
            shading interp
            c = colorbar;
            c.Label.String = 'Surviving Colony Count (a.u.)';
            hold on
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask_quadrat, ...
                'Color', 'yellow', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose_quadrat, ...
                'Color', 'green', 'LineWidth', 1)
            visboundaries(assay.(assay_pattern{i}).img{j}.bw_quadrat, ...
                'Color', 'red', 'LineWidth', 1)
            hold off
            set(gca, 'FontSize', 16)

%             % Plot image
%             h = figure();
%             set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%             imshow(assay.(assay_pattern{i}).img{j}.img_dose, ...
%                 [0 dose_assay(j)], 'colormap', jet(4096))
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
% %             visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose, ...
% %                 'Color', 'green', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_seg, ...
%                 'Color', 'red', 'LineWidth', 0.2)
%             %             visboundaries(temp_bw_dose_peak, 'Color', 'm', 'LineWidth', 1)
%             %             visboundaries(temp_bw_dose_valley, 'Color', 'green', 'LineWidth', 1)
%             %             plot(assay.(assay_pattern{i}).img{j}.y_centroids, ...
%             %                 assay.(assay_pattern{i}).img{j}.x_centroids, 'rx')
%             hold off
%             set(gca, 'FontSize', 16)
%             %             saveas(h, fullfile(path, [fn_assay_img '_imgdose.png']))
% % 
%             h = figure();
%             set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%             imshow(assay.(assay_pattern{i}).img{j}.img_peak, ...
%                 [0 max(assay.(assay_pattern{i}).img{j}.img_peak(:))], ...
%                 'colormap', jet(4096))
%             shading interp
%             c = colorbar;
%             c.Label.String = 'Dose-weighted distance to peak (Gy cm)';
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
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_seg, ...
%                 'Color', 'red', 'LineWidth', 0.2)
%             %             visboundaries(assay.(assay_pattern{i}).img{j}.bw_dose, ...
%             %                 'Color', 'green', 'LineWidth', 1)
%             hold off
%             set(gca, 'FontSize', 16)
%             %             saveas(h, fullfile(path, [fn_assay_img '_imgpeakdist.png']))

%             h = figure();
%             set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%             imshow(img_assay)
%             hold on
%             for row = 1:dxdy_px(1):size(assay.(assay_pattern{i}).img{j}.img_dose,1)
%                 line([1, 2209], [row, row], 'Color', [.7 .7 .7])
%             end
%             for col = 1:dxdy_px(2):size(assay.(assay_pattern{i}).img{j}.img_dose,2)
%                 line([col, col], [1, 3008], 'Color', [.7 .7 .7])
%             end
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_flask, ...
%                 'Color', 'yellow', 'LineWidth', 1)
%             visboundaries(assay.(assay_pattern{i}).img{j}.bw_seg, ...
%                 'Color', 'red', 'LineWidth', 1)
%             hold off
%             set(gca, 'FontSize', 16)
%             %         saveas(h, fullfile(path, [fn_assay_img '_imgassay.tiff']))

        end

    end
end

%% Loop over controls first

N_colonies_ctrl         = [];
PE_ctrl                 = [];
N_colonies_ctrl_quadrat = [];

for i = 1:length(assay.Control.img)

    N_colonies_ctrl = [N_colonies_ctrl assay.Control.img{i}.N_colonies];
    PE_ctrl         = [PE_ctrl assay.Control.img{i}.PE];
    N_colonies_ctrl_quadrat = [N_colonies_ctrl_quadrat; ...
        assay.Control.img{i}.PV.count_quadrat];

end

% [min(PE_ctrl) max(PE_ctrl)].*100
N_colonies_ctrl         = mean(N_colonies_ctrl);
PE_ctrl                 = mean(PE_ctrl);
N_colonies_ctrl_quadrat = mean(N_colonies_ctrl_quadrat);

%% K-fold cross-validation for Poisson regression of quadrat data (training)

close all

lgdstr  = {'Control', 'Dots', 'Dots - peak', 'Dots - valley', ...
    'Stripes', 'Stripes - peak', 'Stripes - valley', 'Open', ...
    'LQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2)})', ...
    'MLQ fit (\it{\lambda_0exp(-\alphaD-\betaD^2+\etaR)})', 'MLQ 95 % CI'};
color   = {'black', 'red', 'green', 'blue'};
N = 1;
K = 4; % number of folds

coeffvariables_LQ_array  = []; % zeros(3, 4, N);
coeffvariables_MLQ_array = []; %zeros(4, 4, N);
mdlcriterion_array       = []; %zeros(2, 2, N);
rmse_array               = []; %zeros(2, 1, N);
mae_array                = []; %zeros(2, 1, N);

for n = 1:N

    % Generate K-fold indices
    cv_indices      = crossvalind('Kfold', 12, K);
    cv_indices_ctrl = crossvalind('Kfold', 4, K);

    assay_test          = {};
    assay_train         = {};
    coeffvariables_LQ   = []; %zeros(3, 4, K);
    coeffvariables_MLQ  = []; %zeros(4, 4, K);
    mdlcriterion        = []; %zeros(2, 2, K);
    mse_values          = []; %zeros(2, 1, K);
    mae_values          = []; %zeros(2, 1, K);

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
        X_D_train    = [];     X_D_test    = [];
        X_PV_R_train = [];     X_PV_R_test = [];
        X_VP_R_train = [];     X_VP_R_test = [];
        X_A_train    = [];     X_A_test    = [];
        X_DP_train   = [];     X_DP_test   = [];
        X_DV_train   = [];     X_DV_test   = [];
        Y_CC_train   = [];     Y_CC_test   = [];
        Y_SF_train   = [];     Y_SF_test   = [];

        % Structure quadrat dataset
        for j = 1:length(assay_pattern)

            dose_vec_train.(assay_pattern{j})               = [];
            peaktovalley_dist_vec_train.(assay_pattern{j})  = [];
            valleytopeak_dist_vec_train.(assay_pattern{j})  = [];
            irradarea_vec_train.(assay_pattern{j})          = [];
            peakdose_vec_train.(assay_pattern{j})           = [];
            valleydose_vec_train.(assay_pattern{j})         = [];
            count_vec_train.(assay_pattern{j})              = [];
            SF_vec_train.(assay_pattern{j})                 = [];

            dose_vec_test.(assay_pattern{j})                = [];
            peaktovalley_dist_vec_test.(assay_pattern{j})   = [];
            valleytopeak_dist_vec_test.(assay_pattern{j})   = [];
            irradarea_vec_test.(assay_pattern{j})           = [];
            peakdose_vec_test.(assay_pattern{j})            = [];
            valleydose_vec_test.(assay_pattern{j})          = [];
            count_vec_test.(assay_pattern{j})               = [];
            SF_vec_test.(assay_pattern{j})                  = [];

            for i = 1:length(assay_train.(assay_pattern{j}).img)

                dose_quadrat_PV              = assay_train.(assay_pattern{j}).img{i}.PV.dose_quadrat;
                peaktovalley_dist_quadrat_PV = assay_train.(assay_pattern{j}).img{i}.PV.peaktovalley_dist_quadrat;
                valleytopeak_dist_quadrat_PV = assay_train.(assay_pattern{j}).img{i}.PV.valleytopeak_dist_quadrat;
                irradarea_quadrat_PV         = assay_train.(assay_pattern{j}).img{i}.PV.irradarea_quadrat;
                peakdose_quadrat_PV          = assay_train.(assay_pattern{j}).img{i}.PV.peakdose_quadrat;
                valleydose_quadrat_PV        = assay_train.(assay_pattern{j}).img{i}.PV.valleydose_quadrat;
                count_quadrat_PV             = assay_train.(assay_pattern{j}).img{i}.PV.count_quadrat;
                SF_quadrat_PV                = count_quadrat_PV ./ N_colonies_ctrl_quadrat;

                dose_vec_train.(assay_pattern{j})              = [dose_vec_train.(assay_pattern{j}) dose_quadrat_PV'];
                peaktovalley_dist_vec_train.(assay_pattern{j}) = [peaktovalley_dist_vec_train.(assay_pattern{j}) peaktovalley_dist_quadrat_PV'];
                valleytopeak_dist_vec_train.(assay_pattern{j}) = [valleytopeak_dist_vec_train.(assay_pattern{j}) valleytopeak_dist_quadrat_PV'];
                irradarea_vec_train.(assay_pattern{j})         = [irradarea_vec_train.(assay_pattern{j}) irradarea_quadrat_PV'];
                peakdose_vec_train.(assay_pattern{j})          = [peakdose_vec_train.(assay_pattern{j}) peakdose_quadrat_PV'];
                valleydose_vec_train.(assay_pattern{j})        = [valleydose_vec_train.(assay_pattern{j}) valleydose_quadrat_PV'];
                count_vec_train.(assay_pattern{j})             = [count_vec_train.(assay_pattern{j}) count_quadrat_PV'];
                SF_vec_train.(assay_pattern{j})                = [SF_vec_train.(assay_pattern{j}) SF_quadrat_PV'];

            end

            for i = 1:length(assay_test.(assay_pattern{j}).img)

                dose_quadrat_PV              = assay_test.(assay_pattern{j}).img{i}.PV.dose_quadrat;
                peaktovalley_dist_quadrat_PV = assay_test.(assay_pattern{j}).img{i}.PV.peaktovalley_quadrat;
                valleytopeak_dist_quadrat_PV = assay_test.(assay_pattern{j}).img{i}.PV.valleytopeak_quadrat;
                irradarea_quadrat_PV         = assay_test.(assay_pattern{j}).img{i}.PV.irradarea_quadrat;
                peakdose_quadrat_PV          = assay_test.(assay_pattern{j}).img{i}.PV.peakdose_quadrat;
                valleydose_quadrat_PV        = assay_test.(assay_pattern{j}).img{i}.PV.valleydose_quadrat;
                count_quadrat_PV             = assay_test.(assay_pattern{j}).img{i}.PV.count_quadrat;
                SF_quadrat_PV                = count_quadrat_PV ./ N_colonies_ctrl_quadrat;

                dose_vec_test.(assay_pattern{j})              = [dose_vec_test.(assay_pattern{j}) dose_quadrat_PV'];
                peaktovalley_dist_vec_test.(assay_pattern{j}) = [peaktovalley_dist_vec_test.(assay_pattern{j}) peaktovalley_dist_quadrat_PV'];
                valleytopeak_dist_vec_test.(assay_pattern{j}) = [valleytopeak_dist_vec_test.(assay_pattern{j}) valleytopeak_dist_quadrat_PV'];
                irradarea_vec_test.(assay_pattern{j})         = [irradarea_vec_test.(assay_pattern{j}) irradarea_quadrat_PV'];
                peakdose_vec_test.(assay_pattern{j})          = [peakdose_vec_test.(assay_pattern{j}) peakdose_quadrat_PV'];
                valleydose_vec_test.(assay_pattern{j})        = [valleydose_vec_test.(assay_pattern{j}) valleydose_quadrat_PV'];
                count_vec_test.(assay_pattern{j})             = [count_vec_test.(assay_pattern{j}) count_quadrat_PV'];
                SF_vec_test.(assay_pattern{j})                = [SF_vec_test.(assay_pattern{j}) SF_quadrat_PV'];

            end

            % Design covariate matrices X
            X_D_train     = [X_D_train     dose_vec_train.(assay_pattern{j})];
            X_PV_R_train  = [X_PV_R_train  peaktovalley_dist_vec_train.(assay_pattern{j})];
            X_VP_R_train  = [X_VP_R_train  valleytopeak_dist_vec_train.(assay_pattern{j})];
            X_A_train     = [X_A_train     irradarea_vec_train.(assay_pattern{j})];
            X_DP_train    = [X_DP_train    peakdose_vec_train.(assay_pattern{j})];
            X_DV_train    = [X_DV_train    valleydose_vec_train.(assay_pattern{j})];
            Y_CC_train    = [Y_CC_train    count_vec_train.(assay_pattern{j})];
            Y_SF_train    = [Y_SF_train    SF_vec_train.(assay_pattern{j})];

            X_D_test      = [X_D_test      dose_vec_test.(assay_pattern{j})];
            X_PV_R_test   = [X_PV_R_test   peaktovalley_dist_vec_test.(assay_pattern{j})];
            X_VP_R_test   = [X_VP_R_test   valleytopeak_dist_vec_test.(assay_pattern{j})];
            X_A_test      = [X_A_test      irradarea_vec_test.(assay_pattern{j})];
            X_DP_test     = [X_DP_test     peakdose_vec_test.(assay_pattern{j})];
            X_DV_test     = [X_DV_test     valleydose_vec_test.(assay_pattern{j})];
            Y_CC_test     = [Y_CC_test     count_vec_test.(assay_pattern{j})];
            Y_SF_test     = [Y_SF_test     SF_vec_test.(assay_pattern{j})];

        end

        X_Dsq_train  = X_D_train.^2;
        X_Rsq_train  = X_PV_R_train.^2;
        X_DRsq_train = X_D_train .* X_Rsq_train;
        X_Dsq_test   = X_D_test.^2;
        X_Rsq_test   = X_PV_R_test.^2;
        X_DRsq_test  = X_D_test .* X_Rsq_test;

        %     [X_D, index_sort]       = sort(X_D);
        %     X_R                             = X_R(index_sort);
        %     X_A                     = X_A(index_sort);
        %     Y_CC                    = Y_CC(index_sort);
        %     [X_Open, index_sort]    = sort(X_Open);
        %     Y_Open                  = Y_Open(index_sort);

        % Perform Poisson regression
        % alpha = 0.24 +- 0.02 Gy-1, beta = 0.019 +- 0.002 Gy-2
        tbl_LQ = table(X_D_train(:), Y_CC_train(:), ...
            'VariableNames', {'D', 'Y_CC'});
        tbl_MLQ = table(X_D_train(:), X_PV_R_train(:), X_A_train(:), Y_CC_train(:), ...
            'VariableNames', {'D', 'R', 'A', 'Y_CC'});
%         tbl_MLQ = table(X_D_train(:), X_P_train(:), Y_CC_train(:), 'VariableNames', {'D', 'P', 'Y_CC'});

        mdl_LQ = fitglm(tbl_LQ, 'Y_CC ~ D + D^2', ...
            'Link', 'log', 'Distribution', 'poisson', 'Intercept', true);
        mdl_MLQ = fitglm(tbl_MLQ, 'Y_CC ~ D + D^2 + R + R^2', ...
            'Link', 'log', 'Distribution', 'poisson', 'Intercept', true);
%         mdl_MLQ = fitglm(tbl_MLQ, 'Y_CC ~ D + D^2 + P', ...
%             'Link', 'log', 'Distribution', 'poisson', 'Intercept', false);

        % Store model coefficient estimates and model criterion estimates
        coeffvariables_LQ(:,:,k)  = mdl_LQ.Coefficients.Variables;
        coeffvariables_MLQ(:,:,k) = mdl_MLQ.Coefficients.Variables;

        mdlcriterion(:,:,k) = [ ...
            mdl_LQ.ModelCriterion.AIC    mdl_LQ.ModelCriterion.BIC; ...
            mdl_MLQ.ModelCriterion.AIC   mdl_MLQ.ModelCriterion.BIC];

        %     G = mdl_LQ.Deviance - mdl_MLQ_R.Deviance;
        %     p = 1 - chi2cdf(G, 1)

        % Create data points for prediction (training)
        D_train          = linspace(min(X_D_train(:)), max(X_D_train(:)), 101);
        R_train          = linspace(min(X_PV_R_train(:)), max(X_PV_R_train(:)), 101);
        [D_mesh, R_mesh] = meshgrid(D_train, R_train);

        % LQ and MLQ predictions (training)
        [yhat_train_LQ, ci_train_LQ]   = predict(mdl_LQ,  [D_mesh(:)]);
        [yhat_train_MLQ, ci_train_MLQ] = predict(mdl_MLQ, [D_mesh(:), R_mesh(:)]);

        % Plot prediction surface
        plot3Dsurf(X_D_train, X_PV_R_train, Y_CC_train, D_mesh, R_mesh, ...
            reshape((yhat_train_LQ), 101, 101), ...
            reshape((yhat_train_MLQ), 101, 101), ...
            'Distance from peak (cm)', ...
            'MLQ fit $(\exp(\lambda_0 -\alpha D - \beta D^2 - \delta A D_p R))$', ...
            fullfile(path, sprintf('SurfPlot_LQvsMLQ_%s.png', num2str(k))))

        % LQ and MLQ predictions (test)
        [yhat_test_LQ, ci_test_LQ]   = predict(mdl_LQ,  [X_D_test(:)]);
        [yhat_test_MLQ, ci_test_MLQ] = predict(mdl_MLQ, [X_D_test(:), X_PV_R_test(:)]);

        % Compute mean squared error and mean absolute error
        mse_values(1,1,k) = mean((Y_CC_test(:) - yhat_test_LQ).^2);
        mse_values(2,1,k) = mean((Y_CC_test(:) - yhat_test_MLQ).^2);
        mae_values(1,1,k) = mean(abs(Y_CC_test(:) - yhat_test_LQ));
        mae_values(2,1,k) = mean(abs(Y_CC_test(:) - yhat_test_MLQ));

    end

    coeffvariables_LQ_array(:,:,n)  = mean(coeffvariables_LQ, 3);
    coeffvariables_MLQ_array(:,:,n) = mean(coeffvariables_MLQ, 3);
    mdlcriterion_array(:,:,n)       = mean(mdlcriterion, 3);
    rmse_array(:,:,n)               = sqrt(mean(mse_values, 3));
    mae_array(:,:,n)                = mean(mae_values, 3);
    disp(n)

end

% Compute 95% confidence interval
significance = 0.000000000000001; % 0.0001; 0.000001 is too unstable for CI bands

[avg_LQ, ci95min_LQ, ci95max_LQ]       = estimateCI(coeffvariables_LQ_array, significance);
[avg_MLQ, ci95min_MLQ, ci95max_MLQ]    = estimateCI(coeffvariables_MLQ_array, significance);
[avg_crit, ci95min_crit, ci95max_crit] = estimateCI(mdlcriterion_array, significance);
[avg_rmse, ci95min_rmse, ci95max_rmse] = estimateCI(rmse_array, significance);
[avg_mae, ci95min_mae, ci95max_mae]    = estimateCI(mae_array, significance);

%%

for j = 2:3 % 1:length(assay_pattern)

    dose_peak_vec.(assay_pattern{j})       = [];
    dose_valley_vec.(assay_pattern{j})     = [];
    peakdist_peak_vec.(assay_pattern{j})   = [];
    peakdist_valley_vec.(assay_pattern{j}) = [];
    peak_peak_vec.(assay_pattern{j})       = [];
    peak_valley_vec.(assay_pattern{j})     = [];

    for i = 1:length(assay.(assay_pattern{j}).img)

        peakdist_peak_vec.(assay_pattern{j}) = [...
            peakdist_peak_vec.(assay_pattern{j}); ...
            mean(assay.(assay_pattern{j}).img{i}.PV.peakdist_peak_quadrat)];
        peakdist_valley_vec.(assay_pattern{j}) = [...
            peakdist_valley_vec.(assay_pattern{j}); ...
            mean(assay.(assay_pattern{j}).img{i}.PV.peakdist_valley_quadrat)];
        peak_peak_vec.(assay_pattern{j}) = [...
            peak_peak_vec.(assay_pattern{j}); ...
            mean(assay.(assay_pattern{j}).img{i}.PV.peak_peak_quadrat)];
        peak_valley_vec.(assay_pattern{j}) = [...
            peak_valley_vec.(assay_pattern{j}); ...
            mean(assay.(assay_pattern{j}).img{i}.PV.peak_valley_quadrat)];

    end

    for k = 1:length(assay_dose)
        irrad_dose.(assay_dose{k}).(assay_pattern{j}).valley.dist = ...
            mean(peakdist_valley_vec.(assay_pattern{j})(4*k-3:4*k));
        irrad_dose.(assay_dose{k}).(assay_pattern{j}).peak.dist = ...
            mean(peakdist_peak_vec.(assay_pattern{j})(4*k-3:4*k));
        irrad_dose.(assay_dose{k}).(assay_pattern{j}).valley.distpeak = ...
            mean(peak_valley_vec.(assay_pattern{j})(4*k-3:4*k));
        irrad_dose.(assay_dose{k}).(assay_pattern{j}).peak.distpeak = ...
            mean(peak_peak_vec.(assay_pattern{j})(4*k-3:4*k));
    end

end

for i = 1:length(assay_dose)
    for j = 2:3 % 1:length(assay_pattern)
        for k = 1:length(GRID_region)
            D_temp  = irrad_dose.(assay_dose{i}).(assay_pattern{j}).(GRID_region{k}).dose;
            R_temp  = irrad_dose.(assay_dose{i}).(assay_pattern{j}).(GRID_region{k}).distpeak;
            [mu_LQ, sigma_LQ]  = LQmdl(D_temp, avg_LQ(1,1), avg_LQ(2,1), ...
                avg_LQ(3,1), avg_LQ(1,2), avg_LQ(2,2), avg_LQ(3,2));
            [mu_MLQ, sigma_MLQ] = MLQmdl(D_temp, R_temp, avg_MLQ(1,1), ...
                avg_MLQ(2,1), avg_MLQ(4,1), avg_MLQ(3,1), avg_MLQ(1,2), ...
                avg_MLQ(2,2), avg_MLQ(4,2), avg_MLQ(3,2));
            disp('-----------------------------------------------------')
            disp([assay_dose{i} ' ' assay_pattern{j} ' ' GRID_region{k}]);
            disp(['LQ: ' num2str(mu_LQ) ' +- ' num2str(sigma_LQ)])
            disp(['MLQ: ' num2str(mu_MLQ) ' +- ' num2str(sigma_MLQ)])
        end
    end
end

%%

varNames = {'Filename', 'Nomial dose (Gy)', ...
    'Surviving colony count (Experimental)', ...
    'Surviving colony count (LQ)', 'Surviving colony count (MLQ)'};

for j = 1:length(assay_pattern)

    filename_vec        = [];
    dose_nominal_vec    = [];
    N_colonies_exp_vec  = [];
    N_colonies_LQ_vec   = [];
    N_colonies_MLQ_vec  = [];

    for i = 1:length(assay.(assay_pattern{j}).img)
        filename_vec = [filename_vec; ...
            string(assay.(assay_pattern{j}).img{i}.filename)];
        dose_nominal_vec = [dose_nominal_vec; ...
            assay.(assay_pattern{j}).img{i}.dose_nom];

        % Observetional count
        N_colonies_exp_vec = [N_colonies_exp_vec; ...
            assay.(assay_pattern{j}).img{i}.N_colonies];

        % Prediction count by LQ and MLQ
        D  = assay.(assay_pattern{j}).img{i}.img_dose_quadrat;
        R  = assay.(assay_pattern{j}).img{i}.img_peakdist_quadrat;
        bw = assay.(assay_pattern{j}).img{i}.bw_quadrat;
        if strcmp((assay_pattern{j}), 'Control')
            yhat_LQ  = exp(avg_LQ(1,1)) .* (ones(size(bw)) .* bw);
            yhat_MLQ = exp(avg_MLQ(1,1)) .* (ones(size(bw)) .* bw);
        else
            yhat_LQ  = exp(avg_LQ(1,1) + avg_LQ(2,1).*D + avg_LQ(3,1).*(D.^2));
            yhat_MLQ = exp(avg_MLQ(1,1) + avg_MLQ(2,1).*D + avg_MLQ(4,1).*(D.^2) + avg_MLQ(3,1).*R);
        end
        N_colonies_LQ_vec  = [N_colonies_LQ_vec; round(sum(yhat_LQ(bw), 'all'))];
        N_colonies_MLQ_vec = [N_colonies_MLQ_vec; round(sum(yhat_MLQ(bw), 'all'))];

    end

    assay_pattern{j}
    if strcmp(assay_pattern{j}, 'Control')
        [h_LQ, p_LQ]   = ttest2(N_colonies_exp_vec, N_colonies_LQ_vec, 'Vartype','unequal');
        [h_MLQ, p_MLQ] = ttest2(N_colonies_exp_vec, N_colonies_MLQ_vec, 'Vartype','unequal');
        [p_LQ, h_LQ]   = ranksum(N_colonies_exp_vec, N_colonies_LQ_vec);
        [p_MLQ, h_MLQ] = ranksum(N_colonies_exp_vec, N_colonies_MLQ_vec);
        [p_LQ, p_MLQ; h_LQ h_MLQ]
    else
        for k = 1:3
            [h_LQ, p_LQ]   = ttest2(N_colonies_exp_vec(4*k-3:4*k), N_colonies_LQ_vec(4*k-3:4*k), 'Vartype','unequal');
            [h_MLQ, p_MLQ] = ttest2(N_colonies_exp_vec(4*k-3:4*k), N_colonies_MLQ_vec(4*k-3:4*k), 'Vartype','unequal');
            [p_LQ, h_LQ]   = ranksum(N_colonies_exp_vec(4*k-3:4*k), N_colonies_LQ_vec(4*k-3:4*k));
            [p_MLQ, h_MLQ] = ranksum(N_colonies_exp_vec(4*k-3:4*k), N_colonies_MLQ_vec(4*k-3:4*k));
            [p_LQ, p_MLQ; h_LQ h_MLQ]
        end
    end
    
%     T = table(filename_vec, dose_nominal_vec, N_colonies_exp_vec, ...
%         N_colonies_LQ_vec, N_colonies_MLQ_vec, 'VariableNames', varNames);
%     writetable(T, fullfile(path, 'data.xlsx'), 'Sheet', assay_pattern{j})

end


%% Poisson regression on solely open data

% Design variable matrix X
X_Open  = [];
Y_Open  = [];

% Structure quadrat dataset
for j = 1:length(assay_pattern)

    dose_vec_train.(assay_pattern{j})      = [];
    count_vec_train.(assay_pattern{j})     = [];
    SF_vec.(assay_pattern{j})        = [];

    for i = 1:length(assay.(assay_pattern{j}).img)

        dose_quadrat_PV  = assay.(assay_pattern{j}).img{i}.PV.dose_quadrat;
        count_quadrat_PV = assay.(assay_pattern{j}).img{i}.PV.count_quadrat;
        dose_vec_train.(assay_pattern{j})  = [dose_vec_train.(assay_pattern{j}) dose_quadrat_PV'];
        count_vec_train.(assay_pattern{j}) = [count_vec_train.(assay_pattern{j}) count_quadrat_PV'];

    end

    % Design covariate matrices X
    if any(strcmp(assay_pattern{j}, {'Control', 'Open'}))
        X_Open  = [X_Open   dose_vec_train.(assay_pattern{j})];
        Y_Open  = [Y_Open   count_vec_train.(assay_pattern{j})];
    end

end

tbl_Open    = table(X_Open(:), Y_Open(:), 'VariableNames', {'D_Open', 'Y_CC_Open'});
mdl_LQ_Open = fitglm(tbl_Open, 'Y_CC_Open ~ D_Open + D_Open^2', ...
    'Link', 'log', 'Distribution', 'poisson', 'Intercept', true);

% Create data points for prediction (training)
D = linspace(min(X_Open(:)), max(X_Open(:)), 1001);

% LQ and MLQ predictions (training)
[yhat_LQ, ci_LQ] = predict(mdl_LQ_Open,  [D(:)]);

figure();
plot(X_Open, Y_Open, 'o', D, yhat_LQ, '-r')
xlabel('Dose (Gy)')
ylabel('Surviving Colony Count')
legend(['Dose response quadrat data $(N=' num2str(length(X_Open)) ')$'], ...
    'LQ fit $(\exp(\lambda_0 -\alpha D - \beta D^2))$', ...
    'Interpreter', 'latex', 'Location', 'northeast')
grid on
set(gca, 'FontSize', 16)

[yhat_LQ, ~] = predict(mdl_LQ_Open,  [X_Open(:)]);
mse_value  = mean((Y_Open(:) - yhat_LQ(:)).^2);
rmse_value = sqrt(mse_value);


coeffvariables_LQ_Open = mdl_LQ_Open.Coefficients.Variables;
mdlcriterion_LQ_Open = [mdl_LQ_Open.ModelCriterion.AIC mdl_LQ_Open.ModelCriterion.BIC];

%% Linear-quadratic (LQ) regression

replicates  = length(dose_assay)/length(unique(dose_assay));
D           = zeros(1, replicates);
SF          = ones(1, replicates);

for K = 1:length(dose_assay)

    D_temp      = dose_assay(ind_sort(K));
    SF_temp     = (assay.Open.img{ind_sort(K)}.N_colonies_uncalib * f(D_temp)) ...
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
D_train   = linspace(0, 10, 1001);
SF_fit  = exp(-alpha.*D_train - beta.*D_train.^2);

figure();
semilogy(D, SF, 'o', D_train, SF_fit, '-r')
xlabel('Dose (Gy)')
ylabel('Surviving Fraction')
legend(['Dose response data $(N=' num2str(length(D)) ')$']', ...
    'LQ fit $(\exp(-\alpha D - \beta D^2))$', ...
    'Interpreter', 'latex', 'Location', 'southwest')
grid on
set(gca, 'FontSize', 16)

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
