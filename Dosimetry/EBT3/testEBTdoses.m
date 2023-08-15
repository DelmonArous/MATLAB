clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path_calib = {'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\New Calibration Data', ...
    'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\EBT3test\Series 2 - best', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\EBT3test\Series 3 - best', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Dosimetry\EBT3 scan Olga\EBT3test\Black'};

%% Variables

% Scanning resolution (in dpi), bit depth and ROI size (in mm)
dpi = 300;
bit = 48;
ROI_size_mm = [4 4];

% EBT3 film dose (in Gy)
dose = [10.0 10.0 10.0 10.0 5.0 5.0 5.0 5.0 7.0 7.0 7.0 7.0];
% dose = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% dose = [10.0 10.0 10.0 10.0];
%dose = [5.0 5.0 5.0 5.0];
% dose = [7.0 7.0 7.0 7.0];
% dose = [0 0 0 0];

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = [1  580];
y_range_px = [35 613];

% Plot spesifications
channel     = {'red', 'green', 'blue', 'gray'};
markerspec  = {'ro', 'gx', 'bs', 'kd'};
linespec    = {'r-', 'g-', 'b-', 'k-'};
lgdstr      = {'Red channel', 'Green channel', 'Blue channel', ...
    'Grayscale image'};

%%
EBT3_calib = getEBT3struct(path_calib{3}, dpi, bit, ROI_size_mm, ...
    x_range_px, y_range_px);

struct_calib.red.PV_control     = 4.5735e+04; % 4.5538e+04; % 
struct_calib.green.PV_control   = 4.2917e+04; % 4.2688e+04; % 
struct_calib.blue.PV_control    = 2.9562e+04; % 2.9367e+04; % 
struct_calib.gray.PV_control    = 4.2237e+04; % 4.2026e+04; %

struct_calib.red.PV_bckg        = 848.7026; % 863.1957;
struct_calib.green.PV_bckg      = 883.3738; % 845.4016;
struct_calib.blue.PV_bckg       = 894.0910; % 822.28;
struct_calib.gray.PV_bckg       = 892.6777; % 847.9192;

% sum = 0;
% for j = 1:length(channel)
% for i = 1:length(EBT3_calib)
%    
%     sum = sum + EBT3_calib{i}.(channel{j}).PV;
%     
% end 
%     channel{j}
% 
% sum = sum / length(EBT3_calib)
% 
% end

figure();
hold on
for j = 1:length(channel)
    
    ODvec       = [];
    netODvec    = [];
    PVvec       = [];
    
    for i = 1:length(EBT3_calib) % length(EBT3_calib)-2
        
        PVvec = [PVvec EBT3_calib{i}.(channel{j}).PV];
        ODvec = [ODvec EBT3_calib{i}.(channel{j}).OD];
        netOD = log10(struct_calib.(channel{j}).PV_control - ...
            struct_calib.(channel{j}).PV_bckg) / ...
            (EBT3_calib{i}.(channel{j}).PV - ...
            struct_calib.(channel{j}).PV_bckg);
        netODvec = [netODvec netOD];
        
    end
    
    % SD and MAD for films
%     [channel{j} ', 5 Gy: SD=' num2str(std(ODvec(5:8))) ', MAD=' num2str(mad(ODvec(5:8)))]
%     [channel{j} ', 7 Gy: SD=' num2str(std(ODvec(9:12))) ', MAD=' num2str(mad(ODvec(9:12)))]
%     [channel{j} ', 10 Gy: SD=' num2str(std(ODvec(1:4))) ', MAD=' num2str(mad(ODvec(1:4)))]
    
%     [channel{j} ', 5 Gy: PV=' num2str((PVvec(5:8)))]
%     [channel{j} ', 7 Gy: PV=' num2str((PVvec(9:12)))]
%     [channel{j} ', 10 Gy: PV=' num2str((PVvec(1:4)))]
    
    % SD and MAD for black
%     [channel{j} ': SD=' num2str(std(ODvec)) ', MAD=' num2str(mad(ODvec))]
%     [channel{j} ', PV=' num2str((PVvec))]
    
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


%%
% figure();
% hold on
% for k = 1:1000
%     
%     ind = randi([1 16],4,1)';
%     
%     for j = 1:length(channel)
%         
%         ODvec       = [];
%         netODvec    = [];
%         MADvec      = [];
%         SDvec       = [];
%         
%         for i = ind % :length(EBT3_calib) % length(EBT3_calib)-2
%             
%             ODvec = [ODvec EBT3_calib{i}.(channel{j}).OD];
%             netOD = log10(struct_calib.(channel{j}).PV_control - ...
%                 struct_calib.(channel{j}).PV_bckg) / ...
%                 (EBT3_calib{i}.(channel{j}).PV - ...
%                 struct_calib.(channel{j}).PV_bckg);
%             netODvec = [netODvec netOD];
%             
%         end
%         
% %         [channel{j} ' ' num2str(std(ODvec)) ' ' num2str(mad(ODvec))]
%         
%         h(j) = plot(k, mad(ODvec), markerspec{j});
%         
%     end
%     
% end
% xlabel('Iteration')
% ylabel('\it{MAD}')
% % xlim([0 101])
% % ylim([11 12])
% legend([h(1), h(2), h(3), h(4)], ...
%     'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
%     'Data, grayscale image', 'Location', 'NorthWest')
% grid on
% set(gca, 'FontSize', 16)
% hold off



% for j = 1:length(channels)
%
%     % Initialize sums and counters
%     sum1_bckg       = 0;    sum2_bckg       = 0;
%     N_bckg    = 0;
%
%     for i = [13 14]
%
%         sum1_bckg = sum1_bckg + ...
%             (EBT3_calib{i}.(channel{j}).PV / ...
%             EBT3_calib{i}.(channel{j}).sigma^2);
%         sum2_bckg = sum2_bckg + ...
%             (1/EBT3_calib{i}.(channel{j}).sigma^2);
%         N_bckg = N_bckg + 1;
%
%     end
%
%     % Calculate mean values
%     struct_calib.(channel{j}).PV_bckg    = sum1_bckg/sum2_bckg;
%
%     % Calculate corresponding standard deviations
%     struct_calib.(channel{j}).sigma_bckg    = sqrt(N_bckg/sum2_bckg);
%
% end
%
% %% Store raw data (dose and netOD) from respective channels into vectors
%
% % Loop over all EBT3 films
% for j = 1:length(channel)
%
%     % Initialize vectors
%     struct_calib.(channel{j}).netOD_vec = [];
%     struct_calib.(channel{j}).sigma_netOD_vec = [];
%
%    for i = 1:(length(EBT3_calib)-N_bckg) % (length(EBT3_calib)-N_bckg) % (N_control+1)
%
%        if i == ...  %  i == 9 || i == 10 % i == 1 || i == 2
%            EBT3_calib{i}.(channel{j}).netOD         = 0;
%            EBT3_calib{i}.(channel{j}).sigma_netOD   = 0;
%        else
%            % Calculate netOD and its corresponding standard deviation
%            EBT3_calib{i}.(channel{j}).netOD = calculateNetOD( ...
%                struct_calib.(channel{j}).PV_control, ...
%                EBT3_calib{i}.(channel{j}).PV, ...
%                struct_calib.(channel{j}).PV_bckg);
%
%            EBT3_calib{i}.(channel{j}).sigma_netOD = calculateSigmaNetOD( ...
%                struct_calib.(channel{j}).PV_control, ...
%                EBT3_calib{i}.(channel{j}).PV, ...
%                struct_calib.(channel{j}).PV_bckg, ...
%                struct_calib.(channel{j}).sigma_control, ...
%                EBT3_calib{i}.(channel{j}).sigma, ...
%                struct_calib.(channel{j}).sigma_bckg);
%        end
%
% %        [EBT3_calib{i}.(channel{j}).netOD, ...
% %            EBT3_calib{i}.(channel{j}).sigma_netOD] = calculateNetOD( ...
% %            EBT3_calib{i}.(channel{j}), struct_calib.(channel{j}));
%
%        % Store calculations
%        struct_calib.(channel{j}).netOD_vec = [ ...
%            struct_calib.(channel{j}).netOD_vec ...
%            EBT3_calib{i}.(channel{j}).netOD];
%        struct_calib.(channel{j}).sigma_netOD_vec = [ ...
%            struct_calib.(channel{j}).sigma_netOD_vec ...
%            EBT3_calib{i}.(channel{j}).sigma_netOD];
%
%    end
%
% end

