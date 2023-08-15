clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path_calib  = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Calibration_v1';
path_ctrl   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Control_v1';
path_bckg   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Background';
path_open   = 'C:\Users\delmo\Desktop\Images\21032022\Open';
path_GRID   = 'C:\Users\delmo\Desktop\Images\21032022\GRID Stripes';


path_dest = 'C:\Users\delmo\Desktop\GRID';


%% Variables

% Scanning resolution (in dpi), bit depth and ROI size (in mm)
dpi = 300; % [1200 300 1200 300];
bit = 48;
bitperchannel = bit/3;
ROI_size_mm = [4 4];

% EBT3 film dose (in Gy)
dose = repmat([0 0.1 0.2 0.5 10 1 2 5]', 1, 8)';
dose = dose(:)';

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = [1 470]; % [1  570];
y_range_px = [10 500]; % ;

% Plot spesifications
channel             = {'red', 'green', 'blue', 'gray'};
% markerspec          = {'ro', 'gx', 'bs', 'kd'};
markerspec          = {{'ro', 'go', 'bo', 'ko'}, {'rs', 'gs', 'bs', 'ks'}};
% linespec            = {'r-', 'g-', 'b-', 'k-'};
linespec            = {{'r-', 'g-', 'b-', 'k-'}, {'r--', 'g--', 'b--', 'k--'}};
lgdstr_profiles     = {'Red channel', 'Green channel', 'Blue channel', ...
    'Grayscale image'};

%% Read EBT3 films

EBT3_calib      = getEBT3struct(path_calib, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px);

EBT3_ctrl       = getEBT3struct(path_ctrl, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px);

img_bckg        = getEBT3struct(path_bckg, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px);

EBT3_open       = getEBT3struct(path_open, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px);

EBT3_GRID    = getEBT3struct(path_GRID, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px);

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
grid on
set(gca, 'FontSize', 16)
hold off

%% Convert experimental EBT3 films into dose by using the calibration

% clc
% close all

channel_opt = 'red';

pattern             = {'open', 'GRID'};
linespec            = {'-r', '-b'};
lgdstr_profiles     = {'Open, $n=5$', 'GRID, $n=5$'};

center_EBT3     = 300;
x1_EBT3         = center_EBT3 - 75;
x2_EBT3         = center_EBT3 + 75;
px_size         = (2.54/dpi);     % in cm/pixel
significance    = 0.0001; % 0.000001 is too unstable for CI bands

% Valley and peak positions
pos_valley  = [55:70 170:250 350:430 530:545];
pos_peak    = [105:135 285:315 465:495];
pos_open    = pos_valley(1):pos_peak(end);

h_profilevec    = [];
dose_open       = [];
dose_valley     = [];
dose_peak       = [];

profile = {};

% Differential response (high, low)
ind.open = {[2], [1 3:5]};
ind.GRID = {[2 5], [1 3:4]};

for i = 1:length(pattern) % 2
    
    img_dose = [];
    for j = 1:length(ind.(pattern{i}))
        
        temp_img_dose = convertEBT3netODtoDose( ...
            eval(sprintf('EBT3_%s', pattern{i})), ...
            channel_opt, fit.(response{j}), ctrl, bckg, ind.(pattern{i}){j});
        img_dose = cat(3, img_dose, temp_img_dose);
        
    end
    struct.(pattern{i}).img_dose = img_dose;
    
    doseprofiles = [];
    
    figure(100)
    hold on
    for j = 1:size(struct.(pattern{i}).img_dose, 3)
        
        % Dose profiles
        temp_img = struct.(pattern{i}).img_dose(1:end, x1_EBT3:x2_EBT3, j); % 60:540
        temp_img(temp_img <= 0) = 0;
        
        temp_doseprofile    = mean(temp_img, 2);
        position            = (1:length(temp_doseprofile)) .* px_size;
        
        % Save peak and valley dose for each profiles of a pattern
        if strcmp('open', pattern{i})
            dose_open = [dose_open mean(temp_doseprofile(pos_open))];
        else
            dose_valley = [dose_valley mean(temp_doseprofile(pos_valley))];
            dose_peak   = [dose_peak mean(temp_doseprofile(pos_peak))];
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
    xlabel('Position (cm)')
    ylabel('Dose (Gy)')
    % xlim([30 length(temp_profile_open)-30])
    %     xlim(([30 length(temp_profile)-30].* px_size))
    xlim([0.5 4.5])
    ylim([0 7])
    legend(h_profilevec, lgdstr_profiles, 'Location', 'NorthEast', ...
        'Interpreter', 'LaTeX')
    grid on
    set(gca, 'FontSize', 16)
    hold off
    
    if strcmp('open', pattern{i})
        [avgDose_open, CI_open] = estimateCI(dose_open, significance);
    else
        [avgDose_valley, CI_valley] = estimateCI(dose_valley, significance);
        [avgDose_peak, CI_peak] = estimateCI(dose_peak, significance);
    end
    
    % Estimate 95% CI bands for the dose profiles
    profile.(pattern{i}).Dose = estimateCIProfileBand( ...
        doseprofiles, position, significance);
    
    %     Normalization
    if i == 1
        norm = mean(profile.(pattern{i}).Dose.avgprofile);
    end
    profile.(pattern{i}).Dose.avgprofile = ...
        profile.(pattern{i}).Dose.avgprofile ./ norm;
    profile.(pattern{i}).Dose.CI95_value = ...
        profile.(pattern{i}).Dose.CI95_value ./ norm;
    
end
