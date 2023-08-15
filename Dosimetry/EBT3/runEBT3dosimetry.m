clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files
path_calib = { ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\New folder\BMF_v2', ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\New folder\BMF_v3'};

path_irrad = { ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\New folder\Radium_v2', ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\New folder\Radium_v3'};

%% Variables

% Scanning resolution (in dpi) and ROI size (in mm)
bit = 48;
dpi = 1200;
ROI_size_mm = [4 4];

% EBT3 film dose (in Gy)
dose = [1.0 1.0 3.0 3.0 5.0 5.0 7.0 7.0 0.0 0.0];

% Plot spesifications
channel     = {'red', 'green', 'blue', 'gray'};
markerspec  = {'ro', 'go', 'bo', 'ko'};
linespec    = {'r-', 'g-', 'b-', 'k-'};
lgdstr      = {'Red channel', 'Green channel', 'Blue channel', ...
    'Grayscale image'};

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = [1   2354];  % [1   2351]
y_range_px = [242 2598];  % [111 2456]

best_channel = {'red', 'green'};
ROI_bb = {{[1260 850], [1170 920], [970 1500], [1080 1370], [1420 940], ...
    [1390 850], [1090 1080], [1100 1020], [1000 1000], [1140 1000], ...
    [1220 1350], [1260 1290], [1230 1550], [1340 1320], [1140 1500]}, ...
    {[1560 1320], [1280 1520], [1160 1600], [1500 1360], [1520 1380], ...
    [1520 1380], [1150 1200], [1150 1200], [1150 1200], [1150 1200], ...
    [1200 1500], [1200 1450], [1210 1500], [1210 1520], [1300 1460], [1280 1440]}};

%%
EBT3_calib = getEBT3struct(path_calib{2}, dpi, bit, ROI_size_mm, ...
    x_range_px, y_range_px);

figure();
hold on
for j = 1:length(channel)
    
    ODvec       = [];
    netODvec    = [];
    
    for i = 1:length(EBT3_calib)-2
        
        ODvec = [ODvec EBT3_calib{i}.(channel{j}).OD];
%         netOD = log10(struct_calib.(channel{j}).PV_control - ...
%             struct_calib.(channel{j}).PV_bckg) / ...
%             (EBT3_calib{i}.(channel{j}).PV - ...
%             struct_calib.(channel{j}).PV_bckg);
%         netODvec = [netODvec netOD];
        
    end
    
    % SD and MAD for films
    [channel{j} ', 1 Gy: SD=' num2str(std(ODvec(1:2))) ', MAD=' num2str(mad(ODvec(1:2)))]
    [channel{j} ', 3 Gy: SD=' num2str(std(ODvec(3:4))) ', MAD=' num2str(mad(ODvec(3:4)))]
    [channel{j} ', 5 Gy: SD=' num2str(std(ODvec(5:6))) ', MAD=' num2str(mad(ODvec(5:6)))]
    [channel{j} ', 7 Gy: SD=' num2str(std(ODvec(7:8))) ', MAD=' num2str(mad(ODvec(7:8)))]
    
    % SD and MAD for black
%     [channel{j} ': SD=' num2str(std(netODvec)) ', MAD=' num2str(mad(netODvec))]
    
    h(j) = plot(dose, ODvec, markerspec{j});
    
end
xlabel('Dose (Gy)')
ylabel('\it{OD}') % 'Transmittance \it{T}'  '\it{netOD}' % '\it{OD}'
xlim([min(dose)-0.5 max(dose)+0.5])
% ylim([11 12])
legend([h(1), h(2), h(3), h(4)], ...
    'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
    'Data, grayscale image', 'Location', 'NorthWest')
grid on
set(gca, 'FontSize', 16)

%% Loop over control and background films to compute mean PV

for k = 2 % 1:length(path_calib)
    
% EBT3_calib = getEBT3struct(path_calib{k}, dpi, bit, ROI_size_mm, ...
%     x_range_px, y_range_px);
struct_calib = {};

for j = 1:length(channel) 
    
    % Initialize sums and counters
    sum1_control    = 0;    sum2_control    = 0;
    sum1_bckg       = 0;    sum2_bckg       = 0;
    N_control = 0;
    N_bckg    = 0;    
    
    for i = [9, 10, 11, 12] % [1, 2, 11, 12]
        
        if i == 9 || i == 10
            
            sum1_control = sum1_control + ...
                (EBT3_calib{i}.(channel{j}).PV / ...
                EBT3_calib{i}.(channel{j}).sigma^2);
            sum2_control = sum2_control + ...
                (1/EBT3_calib{i}.(channel{j}).sigma^2);
            N_control = N_control + 1;
             
        else
            
            sum1_bckg = sum1_bckg + ...
                (EBT3_calib{i}.(channel{j}).PV / ...
                EBT3_calib{i}.(channel{j}).sigma^2);
            sum2_bckg = sum2_bckg + ...
                (1/EBT3_calib{i}.(channel{j}).sigma^2);
            N_bckg = N_bckg + 1;
            
        end 
        
    end
    
    % Calculate mean values
    struct_calib.(channel{j}).PV_control = sum1_control/sum2_control; 
    struct_calib.(channel{j}).PV_bckg    = sum1_bckg/sum2_bckg;
    
    % Calculate corresponding standard deviations
    struct_calib.(channel{j}).sigma_control = sqrt(N_control/sum2_control);
    struct_calib.(channel{j}).sigma_bckg    = sqrt(N_bckg/sum2_bckg);
    
end

%% Store raw data (dose and netOD) from respective channels into vectors 

% Loop over all EBT3 films
for j = 1:length(channel)
    
    % Initialize vectors
    struct_calib.(channel{j}).netOD_vec = [];
    struct_calib.(channel{j}).sigma_netOD_vec = [];
    
   for i = 1:(length(EBT3_calib)-N_bckg) % (length(EBT3_calib)-N_bckg) % (N_control+1)
       
       if i == 9 || i == 10 % i == 1 || i == 2
           EBT3_calib{i}.(channel{j}).netOD         = 0;
           EBT3_calib{i}.(channel{j}).sigma_netOD   = 0;
       else
           % Calculate netOD and its corresponding standard deviation
           EBT3_calib{i}.(channel{j}).netOD = calculateNetOD( ...
               struct_calib.(channel{j}).PV_control, ...
               EBT3_calib{i}.(channel{j}).PV, ...
               struct_calib.(channel{j}).PV_bckg);
           
           EBT3_calib{i}.(channel{j}).sigma_netOD = calculateSigmaNetOD( ...
               struct_calib.(channel{j}).PV_control, ...
               EBT3_calib{i}.(channel{j}).PV, ...
               struct_calib.(channel{j}).PV_bckg, ...
               struct_calib.(channel{j}).sigma_control, ...
               EBT3_calib{i}.(channel{j}).sigma, ...
               struct_calib.(channel{j}).sigma_bckg);
       end
       
%        [EBT3_calib{i}.(channel{j}).netOD, ...
%            EBT3_calib{i}.(channel{j}).sigma_netOD] = calculateNetOD( ...
%            EBT3_calib{i}.(channel{j}), struct_calib.(channel{j}));
       
       % Store calculations
       struct_calib.(channel{j}).netOD_vec = [ ...
           struct_calib.(channel{j}).netOD_vec ...
           EBT3_calib{i}.(channel{j}).netOD];
       struct_calib.(channel{j}).sigma_netOD_vec = [ ...
           struct_calib.(channel{j}).sigma_netOD_vec ...
           EBT3_calib{i}.(channel{j}).sigma_netOD];
 
   end
       
end

%% Linear regression

for j = 1:length(channel)
    
    RsquaredAdjusted_opt = 0;
    
    % Search space for the model parameter n, while keeping it fixed as a
    % constant
    for n = 0:0.1:5
        
        % Design matrix
        X = [struct_calib.(channel{j}).netOD_vec' ...
            (struct_calib.(channel{j}).netOD_vec').^n];
        
        % Fit linear model
        mdl = fitlm(X, dose(1:(length(EBT3_calib)-N_bckg))', ...
            'Intercept', false); % (N_control+1)
        
        % Store fitted model with best adjusted R-squared
        if mdl.Rsquared.Adjusted > RsquaredAdjusted_opt
            RsquaredAdjusted_opt = mdl.Rsquared.Adjusted;
            struct_calib.(channel{j}).mdl = mdl;
            struct_calib.(channel{j}).n = n;
            struct_calib.(channel{j}).RsquaredAdjusted = mdl.Rsquared.Adjusted;
        end
       
    end
    
end

%% Plot

figure();
hold on
for j = 1:length(channel)
    
    % Interpolate model fit to a finer grid    
    struct_calib.(channel{j}).netOD_interp = linspace( ...
        min(struct_calib.(channel{j}).netOD_vec), ...
        max(struct_calib.(channel{j}).netOD_vec), 10000);
    struct_calib.(channel{j}).dose_interp = model1( ...
        struct_calib.(channel{j}).netOD_interp, ...
        struct_calib.(channel{j}).mdl.Coefficients.Estimate(1), ...
        struct_calib.(channel{j}).mdl.Coefficients.Estimate(2), ...
        struct_calib.(channel{j}).n);
    
%     ind = find(struct_calib.(channel{j}).dose_interp > max(dose));
%     struct_calib.(channel{j}).dose_interp( ...
%         struct_calib.(channel{j}).dose_interp > max(dose))  = [];
%     struct_calib.(channel{j}).netOD_interp(ind)             = [];
%     struct_calib.(channel{j}).dose_interp = interp1( ...
%         struct_calib.(channel{j}).netOD_vec, ...
%         struct_calib.(channel{j}).mdl.Fitted, ...
%         struct_calib.(channel{j}).netOD_interp, 'spline', 0);
%     h_fit(j) = plot(struct_calib.(channel{j}).mdl.Fitted,  ...
%         struct_calib.(channel{j}).netOD_vec, linespec{j});    

    % Plot
    h_data(j) = plot(dose(1:(length(EBT3_calib)-N_bckg)), ...
        struct_calib.(channel{j}).netOD_vec, markerspec{j}); % (N_control+1)
    h_fit(j) = plot(struct_calib.(channel{j}).dose_interp,  ...
        struct_calib.(channel{j}).netOD_interp, linespec{j});
     
end
xlabel('Dose (Gy)')
ylabel('\it{netOD}')
xlim([min(dose)-0.5 max(dose)+0.5])
legend([h_data(1), h_data(2), h_data(3), h_data(4), h_fit(4)], ...
    'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
    'Data, grayscale image', 'Fit (\it{a\cdotnetOD+b\cdotnetOD^n})', ...
    'Location', 'NorthWest')
grid on
set(gca, 'FontSize', 14)

%% Convert experimental EBT3 films into dose by using the calibration 

% clc
% close all
EBT3_exp{k} = convertEBT3filmsToDose(path_irrad{k}, ...
    struct_calib.(best_channel{k}), best_channel{k}, ROI_bb{k});


end
