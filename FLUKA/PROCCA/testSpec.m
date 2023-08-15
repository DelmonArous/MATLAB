clear all;
% close all;
fclose('all');
clc;

%% File directories
% pathDest    = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp';
pathMC = { ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu\Dose3D_Xray_180kV_01mmCu.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.3 mm Cu\Dose3D_Xray_180kV_03mmCu.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\200 kV 0.1 mm Cu\Dose3D_Xray_200kV_01mmCu.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\200 kV 0.3 mm Cu\Dose3D_Xray_200kV_03mmCu.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\220 kV 0.1 mm Cu\Dose3D_Xray_220kV_01mmCu.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\220 kV 0.3 mm Cu\Dose3D_Xray_220kV_03mmCu.dat'};

% pathMC = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu - RCC\Dose3D_Xray180kV_01mmCu.dat', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu - RPP\Dose3D_Xray180kV_01mmCu.dat'};

energy          = [180 180 200 200 220 220];    % in kV
thickness_Cu    = [0.1 0.3 0.1 0.3 0.1 0.3];    % in mm

% energy          = [180 180];    % in kV
% thickness_Cu    = [0.1 0.1];    % in mm

%% Coordinate ranges (in cm) for x, y and z
% x = [-25.00 25.00]; % [0 50.00], deltax = -25.00
% y = [-28.30 21.70]; % [0 50.00], deltay = -28.30
% z = [-25.00 2.120]; % [0 27.12], deltaz = -25.00

x = [0 6.0];
y = [0 3.0];
z = [0 3.0];

%% Number of bins in x, y and z dimension, respectively
nbins.x = 113;
nbins.y = 512;
nbins.z = 512;

% width = [(y(2)-y(1))/nbins.y (z(2)-z(1))/nbins.z]; % (x(2)-x(1))/nbins.x
% start = [0 0];

for k = 1:length(pathMC)
    
    %% Extract simulated FLUKA MC dose distribution
    
    fileID = fopen(pathMC{k}, 'r');
    
    formatSpec = '%f %f %f %f %f %f %f %f %f %f';
    sizeA = [10 Inf];
    
    A = fscanf(fileID, formatSpec, sizeA);
    A = A.';
    
    fclose(fileID);
    
    % Reshape Fortran matrix A into a 3-D dose matrix
    B = A.';
    doseMC = reshape(B(1:end-8), [nbins.x nbins.y nbins.z]);
    doseMC = doseMC .* 1.602176462*10^(-7) * 10^9; % in nGy
       
    %% Create a binary mask for target cross section
    [n_cols, n_rows] = meshgrid(1:size(doseMC,2), ...
        1:size(doseMC,3));
    centerX = 256; centerY = 256;
    radius = 256;
    
    bw_o = (n_rows - centerY).^2 + (n_cols - centerX).^2 <= radius.^2;
    bw_o(:,1:256) = 0;
    bw_i = (n_rows - centerY).^2 + (n_cols - centerX).^2 <= radius.^2;
    bw_i(:,256:end) = 0;
    
    bw3D_exit = zeros(size(doseMC,1), size(doseMC,2), size(doseMC,3));
    bw3D_entrance = zeros(size(doseMC,1), size(doseMC,2), size(doseMC,3));
    for slice = 1:nbins.x
        bw3D_exit(slice,:,:) = bw_o;
        bw3D_entrance(slice,:,:) = bw_i;
    end
    
    %%  
    dosedata_o = reshape(doseMC.*bw3D_exit, 1, []);
    dosedata_o(dosedata_o == 0) = []; % < 0.6*10^(-5)
    dosedata_i = reshape(doseMC.*bw3D_entrance, 1, []);
    dosedata_i(dosedata_i == 0) = []; % < 0.6*10^(-5)
    
    %% Compute 95% CI
    SEM = std(dosedata_i)/sqrt(length(dosedata_i));  % Standard Error
    ts = tinv([0.025  0.975], length(dosedata_i)-1); % T-Score
    CI_i{k} = mean(dosedata_i) + ts*SEM;             % Confidence Intervals
    avgdose_i(k) = mean(dosedata_i);
    
    SEM = std(dosedata_o)/sqrt(length(dosedata_o));  % Standard Error
    ts = tinv([0.025  0.975], length(dosedata_o)-1); % T-Score
    CI_o{k} = mean(dosedata_o) + ts*SEM;             % Confidence Intervals
    avgdose_o(k) = mean(dosedata_o);
    
end

%% Plot

colorspec = {'red', 'blue', 'red', 'blue', 'red', 'blue'};

figure();
hold on
for k = 1:length(pathMC)
    
    plot_i(k) = plot(energy(k), avgdose_i(k), '-o', ...
        'MarkerSize', 6, 'LineWidth', 1, ...
        'MarkerEdgeColor', [.2 .2 .2], ...
        'MarkerFaceColor', colorspec{k});
    plot_o(k) = plot(energy(k), avgdose_o(k), '-s', ...
        'MarkerSize', 6, 'LineWidth', 1, ...
        'MarkerEdgeColor', [.2 .2 .2], ...
        'MarkerFaceColor', colorspec{k});
    
    err_i(k) = errorbar(energy(k), avgdose_i(k), abs(diff(CI_i{k})), ...
        'Color', 'k');
    err_o(k) = errorbar(energy(k), avgdose_o(k), abs(diff(CI_o{k})), ...
        'Color', 'k');
     
end
legend([plot_i(1), plot_o(1), plot_i(2), plot_o(2), err_i(1)], ...
    '0.1 mm Cu, In', '0.1 mm Cu, Out', ...
    '0.3 mm Cu, In', '0.3 mm Cu, Out', ...
    '95% CI')
% legend([plot_i(1), plot_o(1), plot_i(2), plot_o(2), err_i(1)], ...
%     '0.1 mm Cu, In, Cylinder', '0.1 mm Cu, Out, Cylinder', ...
%     '0.1 mm Cu, In, Box', '0.1 mm Cu, Out, Box', ...
%     '95% CI')
xlim([min(energy)-5 max(energy)+5])    
xlabel('X-ray beam energy (kV)')
ylabel('Dose (nGy per primary)')
grid on
set(gca, 'FontSize', 14)

% % Loop over slices and plot
% for slice = 10:100 % nbins.x % 40:60 in x,
%     
%     doseMC_temp = doseMC(slice,:,:);
%     
%     % Remove the extra dimension
%     doseMC_temp = squeeze(doseMC_temp);
%     %     doseMC_temp = doseMC_temp .* bw_target;
%     
%     %     % Pad the image and calculate the start offsets
%     %     s = max(size(doseMC_temp') .* width);
%     %     offset = ((size(doseMC_temp') .* width) - s) / 2;
%     %
%     %     % Interpolate to center imageA, padding with zeros
%     %     doseMC_temp = interp2(doseMC_temp, (width(1):width(1):s)/width(1) + ...
%     %         ((size(doseMC_temp,2) * width(1)) - s)/(2 * width(1)), ...
%     %         (width(2):width(2):s)'/width(2) + ...
%     %         ((size(doseMC_temp,1) * width(2)) - s)/(2 * width(2)), 'spline', 0);
%     %
%     %     % Create spatial reference object based on the start and width inputs
%     %     RI = imref2d(size(doseMC_temp), [start(1) + offset(1) start(1) + ...
%     %         offset(1) + size(doseMC_temp,2) * width(1)], [-(start(2) - offset(2) + ...
%     %         size(doseMC_temp,1) * width(2)) -start(2) + offset(2)]);
%     
%     
%     h = figure();
%     % set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%     %     imshow(ind2rgb(gray2ind((imgBG) / 2048, 256), colormap('gray')), RI)
%     %     hold on
%     imshow(doseMC_temp, ...
%         'DisplayRange', [0 10^(floor(log10(max(max(doseMC_temp)))))], ...
%         'ColorMap', colormap(jet(4096)));
%     hold off
%     cb = colorbar;
%     cb.Label.String = 'Dose (nGy per primary)';
%     caxis([0 10^(floor(log10(max(max(doseMC_temp)))))])
%     xlabel('cm')
%     ylabel('cm')
%     set(gca, 'FontSize', 14)
%     
%     %     saveas(h, fullfile(pathDest, sprintf('CTcolorwash_slice%d.png', slice)))
%     
% end
