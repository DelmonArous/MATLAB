clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

% path_calib  = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Calibration_v1';
% path_ctrl   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Control_v1';
% path_bckg   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Background';
% path_open   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Dose Measurement\Open_v1';
% path_GRID   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Dose Measurement\Stripes_v1';
% path_MC     = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v5\X-ray 220 kV spectrum woGRID\Dose maps\DAT', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v5\X-ray 220 kV spectrum wGRID\Dose maps\DAT'};
% path_assay   = 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation_v2';
% path_bwflask = 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks_v2';
%
% path_dest = 'C:\Users\delmo\Desktop\GRID';

path_calib  = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Calibration_v1';
path_ctrl   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Control_v1';
path_bckg   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Background';
path_open   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Dose Measurement\Open_v1';
path_GRID   = 'C:\Users\delmo\Desktop\Jacob\EBT3\310821 - Stripes\Dose Measurement\Stripes_v1';
path_MC     = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v5\X-ray 220 kV spectrum woGRID\Dose maps\DAT', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v5\X-ray 220 kV spectrum wGRID\Dose maps\DAT'};
path_assay   = 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation_v4';
path_bwflask = 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks_v4';

path_dest = 'C:\Users\delmo\Desktop\GRID\Results_v4';


% 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation\17122020';
% 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation\18112019';
% 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation\20112019';

% 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks\17122020'
% 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks\18112019'
% 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks\20112019'

% 'C:\Users\delmo\Desktop\Jacob\Colony Assay A549 Segmentation_v2'
% 'C:\Users\delmo\Desktop\Jacob\Flask Binary Masks_v2'


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
    ROI_size_mm, x_range_px, y_range_px, 0);

EBT3_ctrl       = getEBT3struct(path_ctrl, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);

img_bckg        = getEBT3struct(path_bckg, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);

EBT3_open       = getEBT3struct(path_open, dpi, bitperchannel, ...
    ROI_size_mm, x_range_px, y_range_px, 0);

EBT3_GRID    = getEBT3struct(path_GRID, dpi, bitperchannel, ...
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
lgdstr_profiles     = {'Open, $n=16$', 'GRID, $n=16$'};

center_EBT3     = 253;
x1_EBT3         = center_EBT3 - 100;
x2_EBT3         = center_EBT3 + 100;
px_size         = (2.54/dpi);     % in cm/pixel
significance    = 0.0001; % 0.000001 is too unstable for CI bands

% Valley and peak positions
pos_valley  = [90:130 280:320 475:515];
pos_peak    = [192:212 385:405 567:587];
pos_open    = pos_valley(1):pos_peak(end);

h_profilevec    = [];
dose_open   = [];
dose_valley = [];
dose_peak   = [];

profile = {};

% Differential response (high, low)
ind.open = {[1 4:5 10:14],  [2:3 6:9 15:16]};
ind.GRID = {[1:7 12 14:16], [8:11 13]};

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
        temp_img = struct.(pattern{i}).img_dose(1:end, x1_EBT3:x2_EBT3, j);
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
    xlim([0.25 5.75])
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

%% Plot EBT3 2D dose maps

for i = 1:length(pattern)

    plot2Dmap(struct.(pattern{i}).img_dose, px_size, [0 6], 'Dose (Gy)')

end

%% Read and plot Monte Carlo (MC) simulated 2D dose maps

center_MC      = 200; % size(MC.(sprintf('%s', pattern{j})).MeanDoseMap,1)/2
x1_MC          = center_MC - 150;
x2_MC          = center_MC + 150;

% Valley and peak positions
pos_valley_MC   = [44:108 141:206 239:304 336:373]; % [51:102 149:200 244:297]
pos_peak_MC     = [111:138 209:234 308:331];
pos_open_MC     = pos_valley_MC(1):pos_peak_MC(end);

for i = 1:length(pattern)

    filelist        = getAllFiles(path_MC{i});
    dosemap_sum     = zeros(401);
    flask_count     = 0;

    for j = 1:length(filelist)

        [filepath, filename, ext] = fileparts(filelist{j});
        [DoseMap_MC, RI_MC] = read2DMap(filelist{j});
        DoseMap_MC = DoseMap_MC .* 1.602176462*10^(-7) * 10^9; % in nGy

        % Flip dose maps relative to flask 1
        if contains(filename, "Flask1", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{i})). ...
                (sprintf('Flask%i', j)).Dose = flip(flip(DoseMap_MC, 1), 2);
        elseif contains(filename, "Flask2", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{i})). ...
                (sprintf('Flask%i', j)).Dose = flip(DoseMap_MC, 1);
        elseif contains(filename, "Flask4", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{i})). ...
                (sprintf('Flask%i', j)).Dose = flip(DoseMap_MC, 2);
        else
            MC.(sprintf('%s', pattern{i})). ...
                (sprintf('Flask%i', j)).Dose = DoseMap_MC;
        end

        MC.(sprintf('%s', pattern{i})).(sprintf('Flask%i', j)). ...
            RI = RI_MC;
        MC.(sprintf('%s', pattern{i})).(sprintf('Flask%i', j)). ...
            filename = filename;

        dosemap_sum = dosemap_sum + ...
            MC.(sprintf('%s', pattern{i})).(sprintf('Flask%i', j)).Dose;
        flask_count = flask_count + 1;

        % Plot
        %         h = figure();
        %         set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %         imshow( ...
        %             MC.(sprintf('%s', pattern{j})).(sprintf('Flask%i', i)).Dose, ...
        %             MC.(sprintf('%s', pattern{j})).(sprintf('Flask%i', i)).RI, ...
        %             [], 'colormap', jet(256))
        %         c = colorbar;
        %         c.Label.String = 'Dose (nGy/primary)'; % 'Dose (%)'
        %         xlabel('cm'); ylabel('cm');
        %         set(gca, 'FontSize', 16)
        %         shading interp
        %     saveas(h, fullfile(destpath, sprintf('%s_Flask%i.png', GRID{k}, i)))
        %     writematrix(struct.(sprintf('%s', GRID{k})).(sprintf('Flask%i', i)), ...
        %         fullfile(destpath, sprintf('%s_Flask%i.txt', GRID{k}, i)))

    end

    % Average dose map over all four flasks
    MC.(sprintf('%s', pattern{i})).MeanDoseMap = dosemap_sum./flask_count;
%     MC.(sprintf('%s', pattern{i})).MeanDoseMap = ...
%         imrotate(MC.(sprintf('%s', pattern{i})).MeanDoseMap, 90);
    MC.(sprintf('%s', pattern{i})).RI =  ...
        MC.(sprintf('%s', pattern{i})).(sprintf('Flask%i', 3)).RI;
%     temp_WorldLimits = MC.(sprintf('%s', pattern{i})).RI.XWorldLimits;
%     MC.(sprintf('%s', pattern{i})).RI.XWorldLimits = ...
%         MC.(sprintf('%s', pattern{i})).RI.YWorldLimits;
%     MC.(sprintf('%s', pattern{i})).RI.YWorldLimits = temp_WorldLimits;

    % Plot
    h = figure();
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(MC.(sprintf('%s', pattern{i})).MeanDoseMap, ...
        MC.(sprintf('%s', pattern{i})).RI, [], 'colormap', jet(256))
    c = colorbar;
    c.Label.String = 'Dose (nGy/primary)'; % 'Dose (%)'
    xlabel('cm'); ylabel('cm');
    grid on
    set(gca, 'FontSize', 26)
    shading interp

    % Spatial resolution in cm/pixel
    RI_MC       = MC.(sprintf('%s', pattern{i})).Flask3.RI;
    pxsize_MC   = [ diff(RI_MC.YWorldLimits) / size(MC.(sprintf('%s', pattern{i})).MeanDoseMap,1) ...
        diff(RI_MC.XWorldLimits) / size(MC.(sprintf('%s', pattern{i})).MeanDoseMap,2) ];

    % Dose profile
    temp_doseprofile = MC.(sprintf('%s', pattern{i})).MeanDoseMap(1:end, x1_MC:x2_MC);
    temp_doseprofile(temp_doseprofile <= 0) = 0;
    temp_doseprofile = mean(temp_doseprofile, 2);
    profile.(pattern{i}).Dose.data = temp_doseprofile;
    temp_doseprofile = smoothdata(temp_doseprofile, 'movmean', 3);
    position = (1:length(temp_doseprofile)) .* pxsize_MC(1);

    if i == 1
        norm_dose = mean(temp_doseprofile);
    end
    temp_doseprofile = temp_doseprofile ./ norm_dose;
    profile.(pattern{i}).Dose.data = profile.(pattern{i}).Dose.data ./ norm_dose;

    if i == 2
        temp_doseprofile(pos_valley_MC) = 2.0 .* temp_doseprofile(pos_valley_MC);
        profile.(pattern{i}).Dose.data(pos_valley_MC) = 2.0 .* profile.(pattern{i}).Dose.data(pos_valley_MC);
    end

    if strcmp('open', pattern{i})
        %         dose_open = mean(temp_profile(pos_open));
    else
        dose_valley_MC = mean(temp_doseprofile(pos_valley_MC));
        dose_peak_MC   = mean(temp_doseprofile(pos_peak_MC));
    end

    profile.(pattern{i}).Dose.avgprofile_MC     = temp_doseprofile;
    profile.(pattern{i}).Dose.position_MC       = position;

    dlmwrite(['MCresults_' pattern{i} '.dat'], ...
        [profile.(pattern{i}).Dose.avgprofile_MC, ...
        profile.(pattern{i}).Dose.position_MC']);

    %     fid = fopen(['MCresults_' pattern{i} '.dat'], 'w');
    %     % Write headers
    %     fprintf(fid, 'MC, Position (cm)\n');
    %     % Write data.
    %     fprintf(fid, '%f, %f', ...
    %         [profile.(pattern{i}).Dose.avgprofile_MC, ...
    %         profile.(pattern{i}).Dose.position_MC']);
    %     fclose(fid);


end

%% 95% CI and FLUKA MC for the dose profiles

lgdstr  = {'EBT3 mean profile, Open, $n=16$', 'EBT3 95\% CI, Open', ...
    'FLUKA MC data, Open', 'FLUKA MC mean profile, Open', ...
    'EBT3 mean profile, GRID, $n=16$', 'EBT3 95\% CI, GRID', ...
    'FLUKA MC data, GRID', 'FLUKA MC mean profile, GRID'};
facecolor_CI95  = {[1.0 0.8 0.8], [0.3010 0.7450 0.9330]};
color_CI95      = {'red', 'blue'};
color_MC        = {[0 0 0], [0.5 0.5 0.5]};

h_vec = [];

for i = 1:length(pattern)

    figure(200)
    hold on
    h_95CIband = fill(profile.(pattern{i}).Dose.CI95_pos, ...
        profile.(pattern{i}).Dose.CI95_value, color_CI95{i}, ...
        'FaceColor', facecolor_CI95{i}, 'EdgeColor', 'none');
    h_95CImean = plot(profile.(pattern{i}).Dose.position, ...
        profile.(pattern{i}).Dose.avgprofile, ...
        linespec{i}, 'LineWidth', 1.0);
    h_MCdata = plot(profile.(pattern{i}).Dose.position_MC - 0.27, ...
        profile.(pattern{i}).Dose.data, '.', ...
        'Color', color_MC{i}, 'LineWidth', 1.0);
    h_MCmean = plot(profile.(pattern{i}).Dose.position_MC - 0.27, ...
        profile.(pattern{i}).Dose.avgprofile_MC, '-', ...
        'Color', color_MC{i}, 'LineWidth', 1.0);
    hold off
    h_vec = [h_vec h_95CImean h_95CIband h_MCdata h_MCmean];

end
legend(h_vec, lgdstr, 'Location', 'NorthEast', 'NumColumns', 2, ...
    'Interpreter', 'LaTeX')
xlabel('Position (cm)')
ylabel('Normalized Dose')
xlim([0.25 5.75])
ylim([0 1.4])
% ylim([0 0.3*10^(-5)])
grid on
set(gca, 'FontSize', 19)

%%

px_size_assay   = 2.54/1200;    % in cm/px ; 1200 dpi assay
cells_seeded    = 30000;        % no. of seeded cells
assay_pattern   = {'ctrl',  'GRID', 'open'};
lgdstr  = {'Mean predicted survival', '95\% CI predicted', ...
    'Mean observed survival', '95\% CI observed', 'RPD'};
assay_dose      = {'Dose2Gy', 'Dose5Gy', 'Dose10Gy'};
lgd_title = {};
y_limits_SF  = {};
y_limits_RPE  = {};

facecolor_CI95_assay    = {[0.4 0.4 0.4],   [0.7 0.7 0.7]};
color_CI95_assay        = {[0 0 0],         [0.5 0.5 0.5]};
y_limits_SF.open        = {[0.5 1.25],     [0.15 1.15],  [0.0 0.9]}; % {[0.5 0.75],       [0.15 0.3],  [0.0 0.04]};
y_limits_SF.GRID        = {[0.5 1.25],     [0.15 1.15],  [0.0 0.9]};
y_limits_RPE.open       = {[-20 20],        [-30 70],     [-75 250]}; % {[-10 20],        [-15 30],     [-20 35]};
y_limits_RPE.GRID       = {[-20 20],        [-30 70],     [-75 250]};
dose_scale              = [2/5 1 2];
dose_assay              = repmat([2 5 10 2 5 10]', 1, 4)';
dose_assay              = dose_assay(:)';
dose_open               = unique(dose_assay);
dose_valley             = [0.36 0.9 1.8]; %  dose_open .* 0.18;
dose_peak               = [1.64 4.1 8.2]; % dose_open .* 0.82;
[~, ind_sort]           = sort(dose_assay, 'ascend');

folderlist_assay    = getAllFolders(path_assay);
folderlist_bwflask  = getAllFolders(path_bwflask);

% Calibration function: f = a + b*D + c*D^2
D_fit       = linspace(0, max(dose_open), 1001);
% f           = 1.7258 - 0.3130.*D_fit + 0.0143.*D_fit.^2;
f           = 1.7530 - 0.3530.*D_fit + 0.0182.*D_fit.^2; 

close all

dx_cm           = 0.10;                         % in cm
dx_px           = round(dx_cm/px_size_assay);   % in px
lgd_title.open  = { 'Open irradiation, 2.0 Gy dose', ...
                    'Open irradiation, 5.0 Gy dose', ...
                    'Open irradiation, 10.0 Gy dose'};
lgd_title.GRID  = { 'GRID irradiation, 2.0 Gy nominal dose', ...
                    'GRID irradiation, 5.0 Gy nominal dose', ...
                    'GRID irradiation, 10.0 Gy nominal dose'};
% lgd_title.open  = {['2.0 Gy dose, $\Delta x = $ ' num2str(dx_cm*10) ' mm'], ...
%     ['5.0 Gy dose, $\Delta x = $ ' num2str(dx_cm*10) ' mm'], ...
%     ['10.0 Gy dose, $\Delta x = $ ' num2str(dx_cm*10) ' mm']};
% lgd_title.GRID  = {['1.6 Gy peak dose, $\Delta x = $ ' num2str(dx_cm*10) ' mm'], ...
%     ['4.1 Gy peak dose, $\Delta x = $ ' num2str(dx_cm*10) ' mm'], ...
%     ['8.2 Gy peak dose, $\Delta x = $ ' num2str(dx_cm*10) ' mm']};

for i = 1:length(folderlist_assay)

    filelist_assay      = getAllFiles(folderlist_assay{i});
    filelist_bwflask    = getAllFiles(folderlist_bwflask{i});

    for j = 1:length(filelist_assay)

        [~, name, ext]  = fileparts(filelist_assay{j});
        data            = readtable(filelist_assay{j});
        y_centroid      = sort(round(data.CentroidY_Coordinate_px_));
        bw_flask        = load(filelist_bwflask{j});

        row_weights     = sum(bw_flask, 2);
        row_weights     = row_weights ./ max(row_weights);

        colony_per_pxrow = [];
        for k = 1:size(bw_flask, 1)
            colony_per_pxrow = [colony_per_pxrow; ...
                length(y_centroid(y_centroid == k))];
        end
        colony_per_pxrow = colony_per_pxrow ./ row_weights;
        colony_per_pxrow(isinf(colony_per_pxrow)) = 0;
        colony_per_pxrow(isnan(colony_per_pxrow)) = 0;

        n_bands = ceil(size(bw_flask,1)/dx_px);
        dummy_zeros = zeros(n_bands*dx_px - size(bw_flask, 1), 1);
        colony_per_pxrow    = cat(1, colony_per_pxrow, dummy_zeros);
        colony_count = sum(reshape(colony_per_pxrow, ...
            length(colony_per_pxrow)/n_bands, n_bands));
        position = linspace(0, size(bw_flask, 1)*px_size_assay, n_bands); 

        %         size(colony_per_pxrow)
        %         size(colony_count)
        %         size(position)

        %         lower_vec = [];
        %         upper_vec = [];
        %         count_vec = [];
        %         for k = 1:(size(bw_flask,1)/dx)
        %
        %             if k == 1
        %                 lower = 1;
        %             else
        %                 lower = (k - 1) * dx;
        %             end
        %             upper = k * dx;
        %             logical_index = y_centroid >= lower & y_centroid < upper;
        %
        %             lower_vec = [lower_vec lower];
        %             upper_vec = [upper_vec upper];
        %             count_vec = [count_vec sum(logical_index)];
        %
        %         end

        assay.(assay_pattern{i}).img{j}.filename    = name;
        assay.(assay_pattern{i}).img{j}.N_colonies  = length(y_centroid);
        assay.(assay_pattern{i}).img{j}.PE          = ...
            (length(y_centroid)*f(1)) / cells_seeded;
        assay.(assay_pattern{i}).img{j}.position     = position;
        assay.(assay_pattern{i}).img{j}.colony_count = colony_count;
        assay.(assay_pattern{i}).img{j}.avg_count    = median(colony_count);
        if strcmp(assay_pattern{i}, 'ctrl')
            assay.(assay_pattern{i}).img{j}.dose     = 0;
        else
            assay.(assay_pattern{i}).img{j}.dose     = dose_assay(j);
        end

        % assay.(assay_pattern{i}).img{j}.position   = (lower + upper) ./ 2;
        % assay.(assay_pattern{i}).img{j}.count_vec  = count_vec;
        % assay.(assay_pattern{i}).img{j}.profile_PE = count_vec ./ cells_seeded;

    end

end

%% Loop over controls first

N_colonies_ctrl = [];
avg_count_ctrl  = [];
PE_ctrl         = [];
for j = 1:length(assay.ctrl.img)

    N_colonies_ctrl = [N_colonies_ctrl assay.ctrl.img{j}.N_colonies];
    avg_count_ctrl  = [avg_count_ctrl assay.ctrl.img{j}.avg_count];
    PE_ctrl         = [PE_ctrl assay.ctrl.img{j}.PE];

end
[min(PE_ctrl) max(PE_ctrl)].*100
N_colonies_ctrl = mean(N_colonies_ctrl);
avg_count_ctrl  = median(avg_count_ctrl);
PE_ctrl         = mean(PE_ctrl);

%     %% Loop over open 10 Gy irradiations and replace 'correct' colony count
%     for j = 25:length(assay.open.img)
%
%         assay.open.img{j}.N_colonies = N_colony_10Gy(j-24);
%         assay.open.img{j}.PE         = N_colony_10Gy(j-24)/cells_seeded;
%
%     end

%% Estimate and plot EBT3 2D surviving fraction (SF) maps

replicates  = length(dose_assay)/length(unique(dose_assay));
D           = zeros(1, replicates); % zeros(1, length(assay.ctrl.img));
% D           = [D dose_assay];
SF          = ones(1, length(assay.ctrl.img)); % [1]; % ones(1, length(assay.ctrl.img));            

for k = 1:length(dose_assay) % -1

%     assay.open.img{ind_sort(k)}.filename
    D_temp      = dose_assay(ind_sort(k));
    ind_temp    = find(D_fit == D_temp);
    SF_temp     = (assay.open.img{ind_sort(k)}.N_colonies * f(ind_temp)) ...
        / (N_colonies_ctrl * f(1)) ;  % (cells_seeded * PE_ctrl)]; % 

    assay.open.img{ind_sort(k)}.SF = SF_temp;
    SF      = [SF SF_temp];    
    D       = [D D_temp];

%     D = [D repmat(5*dose_scale(k), 1, replicates)];
    %             SF_sum = 0;
 
%     for j = 1:replicates
% 
%         assay.open.img{(k-1)*replicates+j}.filename
%         SF_temp = assay.open.img{(k-1)*replicates+j}.N_colonies / ...
%             N_colonies_ctrl; % (cells_seeded * PE_ctrl);
%         SF      = [SF SF_temp];
% %         SF_test(k+1) = SF_test(k+1) + SF_temp;
% 
%     end

    %         SF_test(k+1) = SF_test(k+1)/replicates;

end

% LQ regression fit
tbl_LQ  = table(D', log(SF)', 'VariableNames', {'D', 'SF'});
mdl_LQ  = fitlm(tbl_LQ, 'SF ~ D + D^2', 'Intercept', false);

% Store alpha and beta LQ-parameters
alpha   = abs(mdl_LQ.Coefficients.Estimate(1));
beta    = abs(mdl_LQ.Coefficients.Estimate(2));
% alpha       = 0.07; % Bjørg: 0.074 ± 0.027 in Gy-1; Hilde: 0.40 ± 0.05 in Gy-1; 0.07(3) Gy-1 and 0.036(3) Gy-2.
% beta        = 0.0285; % Bjørg: 0.026 ± 0.003 in Gy-2; Hilde: 0.031 ± 0.005 in Gy-2;
SF_fit  = exp(-alpha.*D_fit - beta.*D_fit.^2);
%SF_fit = exp(-0.10*D_fit - 0.0027.*D_fit.^2);

figure();
semilogy(D, SF, 'o', D_fit, SF_fit, '-r') % D, exp(mdl_LQ.Fitted)
xlabel('Dose (Gy)')
ylabel('Surviving fraction')
legend('Dose response data', 'LQ fit (\it{exp(-\alphaD-\beta D^2)})', ...
    'Location', 'southwest')
grid on
set(gca, 'FontSize', 16)

for k = 1:length(dose_open)

    for i = 1:length(pattern)

        % Estimate predicted survival
        struct.(pattern{i}).(assay_dose{k}).img_SF = LQ( ...
            dose_scale(k).*struct.(pattern{i}).img_dose, alpha, beta);

        % Plot predicted survival map
        plot2Dmap(struct.(pattern{i}).(assay_dose{k}).img_SF, ...
            px_size, [0 1], 'Surviving fraction')

    end

end

%%

for k = 1:length(unique(dose_assay))

    h_profilevec    = [];
    SF_open         = [];
    SF_valley       = [];
    SF_peak         = [];

    for i = 1:length(pattern)

        SFprofiles = [];

%         figure(500+k)
%         hold on
        for j = 1:size(struct.(pattern{i}).(assay_dose{k}).img_SF, 3)

            % SF profiles
            temp_img = struct.(pattern{i}).(assay_dose{k}). ...
                img_SF(1:end, x1_EBT3:x2_EBT3, j);
            temp_img(temp_img <= 0) = 0;
            temp_img(temp_img >= 1) = 1;

            temp_SFprofile  = mean(temp_img, 2);
            position        = (1:length(temp_SFprofile)) .* px_size;

            % Save peak and valley SF for each profiles of a pattern
            if strcmp('open', pattern{i})
                SF_open = [SF_open mean(temp_SFprofile(pos_open))];
            else
                SF_valley = [SF_valley mean(temp_SFprofile(pos_valley))];
                SF_peak   = [SF_peak mean(temp_SFprofile(pos_peak))];
            end

            % Store SF profile
            SFprofiles = [SFprofiles; temp_SFprofile'];

            % Plot 1D profile across estimated EBT3 dose maps
            % %                 h = plot(position, temp_SFprofile, linespec{i}, 'LineWidth', 1.0);
            %     yline(dose_valley, '--', 'Valley')
            %     yline(dose_peak, '--', 'Peak')
            %     yline(dose_open, '-.', 'Open')

        end
        % %             h_profilevec = [h_profilevec h];
        % %             xlabel('Position (cm)')
        % %             ylabel('Surviving fraction')
        % %             % xlim([30 length(temp_profile_open)-30])
        % %             %     xlim(([30 length(temp_profile)-30].* px_size))
        % %             xlim([0.25 5.75])
        % %             %         ylim([0.2 1])
        % %             lgd = legend(h_profilevec, lgdstr_profiles, ...
        % %                 'Location', 'SouthEast', 'Interpreter', 'LaTeX');
        % %             title(lgd, lgd_title{k})
        % %             grid on
        % %             set(gca, 'FontSize', 16)
        % %             hold off

        if strcmp('open', pattern{i})
            [avgSF_open, CI_open]       = estimateCI(SF_open, significance);
        else
            [avgSF_valley, CI_valley]   = estimateCI(SF_valley, significance);
            [avgSF_peak, CI_peak]       = estimateCI(SF_peak, significance);
        end

        % Estimate 95% CI bands for the SF profiles
        profile.(pattern{i}).SF_LQ.(assay_dose{k}) = ...
            estimateCIProfileBand(SFprofiles, position, significance);

    end

end

%%

shift = -0.075; % 0.135

for k = 1:length(dose_open)

    for i = 2:length(assay_pattern) %  [3, 2] % length(assay_pattern):-1:2

        profile_count_vec   = [];

        for j = 1:replicates

            assay.(assay_pattern{i}).img{ind_sort((k-1)*replicates+j)}.filename
            position = shift + ...
                assay.(assay_pattern{i}).img{ind_sort((k-1)*replicates+j)}.position;
            profile_count_vec = [profile_count_vec; ...
                assay.(assay_pattern{i}).img{ind_sort((k-1)*replicates+j)}.colony_count];
            %             avg_count_vec = [avg_count_vec; ...
            %                 assay.(assay_pattern{i}).img{j}.avg_count];

        end

        % Estimate 95% CI bands for the SF profiles, per dose
        profile.(assay_pattern{i}).SF_assay.(assay_dose{k}) = ...
            estimateCIProfileBand(profile_count_vec, position, 0.5);

        % Normalization; estimate SF
        if strcmp('open', assay_pattern{i})
            ind_temp    = find(D_fit == dose_open(k));
            profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile = ...
                (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile ...
                .* f(ind_temp)) ./ (avg_count_ctrl .* f(1));
            profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value = ...
                (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value ...
                .* f(ind_temp)) ./ (avg_count_ctrl .* f(1));
        else
            temp_pos = ...
                profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).position;
            temp_pos_CI = ...
                profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_pos;
            ind_valley = (temp_pos > 0.42 & temp_pos < 1.50) | ...
                (temp_pos > 1.99 & temp_pos < 3.09) | ...
                (temp_pos > 3.56 & temp_pos < 4.67) | ...
                (temp_pos > 5.135);
            ind_peak = ~ind_valley;
            ind_valley_CI = (temp_pos_CI > 0.42 & temp_pos_CI < 1.50) | ...
                (temp_pos_CI > 1.99 & temp_pos_CI < 3.09) | ...
                (temp_pos_CI > 3.56 & temp_pos_CI < 4.67) | ...
                (temp_pos_CI > 5.135);
            ind_peak_CI = ~ind_valley_CI;

            ind_temp_valleydose = find(D_fit == dose_valley(k));
            ind_temp_peakdose   = find(D_fit == dose_peak(k));

            % Valley
            profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile(ind_valley) = ...
                (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile(ind_valley) ...
                .* f(ind_temp_valleydose)) ./ (avg_count_ctrl .* f(1));
            profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value(ind_valley_CI) = ...
                (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value(ind_valley_CI) ...
                .* f(ind_temp_valleydose)) ./ (avg_count_ctrl .* f(1));

            % Peak
            profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile(ind_peak) = ...
                (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile(ind_peak) ...
                .* f(ind_temp_peakdose)) ./ (avg_count_ctrl .* f(1));
            profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value(ind_peak_CI) = ...
                (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value(ind_peak_CI) ...
                .* f(ind_temp_peakdose)) ./ (avg_count_ctrl .* f(1));

        end


%         % Normalization; estimate SF
%         profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile = ...
%             (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).avgprofile ...
%             .* f(ind_temp)) ./ (avg_count_ctrl .* f(1));
%         profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value = ...
%             (profile.(assay_pattern{i}).SF_assay.(assay_dose{k}).CI95_value ...
%             .* f(ind_temp)) ./ (avg_count_ctrl .* f(1));

    end

end

%% 95% CI of SF profiles for observed vs. predicted

for k = 1:length(assay_dose)

    for i = 1:length(pattern)

        h_vec = [];

        temp_SF_LQ = interp1( ...
            profile.(pattern{i}).SF_LQ.(assay_dose{k}).position, ...
            profile.(pattern{i}).SF_LQ.(assay_dose{k}).avgprofile, ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).position);

%         RPE = 100 .* abs(temp_SF_LQ - ...
%             profile.(pattern{i}).SF_assay.(assay_dose{k}).avgprofile) ./ ...
%             ((profile.(pattern{i}).SF_assay.(assay_dose{k}).avgprofile + ...
%             temp_SF_LQ)./2);
        RPD = 100 .* (temp_SF_LQ - ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).avgprofile) ./ ...
            ( (profile.(pattern{i}).SF_assay.(assay_dose{k}).avgprofile + ...
            temp_SF_LQ)./2 );

        if strcmp('open', pattern{i})
            temp_pos = ...
                profile.(pattern{i}).SF_assay.(assay_dose{k}).position;
            ind_open = temp_pos > 0.5 & temp_pos < 5.1;
            pattern{i}
            mean(RPD(ind_open))
        else
            temp_pos = ...
                profile.(pattern{i}).SF_assay.(assay_dose{k}).position;
            ind_valley = (temp_pos > 0.55 & temp_pos < 1.35) | ...
                (temp_pos > 2.1 & temp_pos < 2.95) | ...
                (temp_pos > 3.7 & temp_pos < 4.55);
            ind_peak = (temp_pos > 1.55 & temp_pos < 1.95) | ...
                (temp_pos > 3.12 & temp_pos < 3.5) | ...
                (temp_pos > 4.7 & temp_pos < 5.1);
            pattern{i}
            mean(RPD(ind_valley))
            mean(RPD(ind_peak))
        end

        h = figure();
        set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        hold on

        yyaxis left
        % Plot predicted response
        h_95CIband_predicted = fill( ...
            profile.(pattern{i}).SF_LQ.(assay_dose{k}).CI95_pos, ...
            profile.(pattern{i}).SF_LQ.(assay_dose{k}).CI95_value, ...
            color_CI95{i}, 'FaceColor', facecolor_CI95{i}, ...
            'EdgeColor', 'none');
        h_95CImean_predicted = plot( ...
            profile.(pattern{i}).SF_LQ.(assay_dose{k}).position, ...
            profile.(pattern{i}).SF_LQ.(assay_dose{k}).avgprofile, ...
            linespec{i}, 'LineWidth', 1.0);
        % Plot observed response
        h_95CIband_observed = fill( ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).CI95_pos, ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).CI95_value, ...
            color_CI95_assay{i}, 'FaceColor', facecolor_CI95_assay{i}, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5);
        h_95CImean_observed = plot( ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).position, ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).avgprofile, ...
            '-', 'Color', color_CI95_assay{i}, 'LineWidth', 1.0);
        ylabel('Surviving fraction')
        ylim(y_limits_SF.(pattern{i}){k})
        set(gca, 'FontSize', 26) %, 'YScale', 'log')

        % Plot relative percentage difference
        yyaxis right
        h_RPE = plot( ...
            profile.(pattern{i}).SF_assay.(assay_dose{k}).position, ... 
            RPD, '--k');
        ylabel('Relative percentage deviation (%)')
        ylim(y_limits_RPE.(pattern{i}){k})
        yline(0, '-', 'RPD = 0', 'LineWidth', 1, 'FontSize', 18);
        hold off
        h_vec = [h_vec h_95CImean_predicted h_95CIband_predicted ...
            h_95CImean_observed h_95CIband_observed h_RPE];

        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';

        xlabel('Position (cm)')
        xlim([0.25 5.75])
        lgd = legend(h_vec, lgdstr, 'Location', 'NorthEast', ...
            'NumColumns', 3, 'Interpreter', 'LaTeX');
        title(lgd, lgd_title.(pattern{i}){k})
%         grid on
        set(gca, 'FontSize', 26)
        hold off

        %             saveas(h, fullfile(path_dest, sprintf( ...
        %                 'SFprofiles_PredictedVsObserved_%s_%s_%smm.png', ...
        %                 pattern{i}, assay_dose{k}, num2str(dx_cm*100))))


    end

end

