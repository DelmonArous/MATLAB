function [LET_d, LET_f] = LETaverageEstimate(sourcepath, n_detectors, depth)

%% Store all LET data files
filelist = getAllFiles(sourcepath)

%% Read data

% Initialization
LETbin_start = []; LETbin_end = []; LETbin_width = []; LETpoint = [];
value = [];
f = []; d = []; LET_f = []; LET_d = [];

formatSpec = '%f %f %f %f';
dataarray_size = [4 Inf];

% Import data into arrays
for i = 1:n_detectors
    
    fileID = fopen(filelist{i}, 'r');
    fgetl(fileID); fgetl(fileID); % skip the first two lines
    dataarray = fscanf(fileID, formatSpec, dataarray_size);
    dataarray =  dataarray';
    
    LETbin_start(:,i) = dataarray(:,1);
    LETbin_end(:,i) = dataarray(:,2);
    value(:,i) = dataarray(:,3);
    
    LETbin_width(:,i) = abs(LETbin_end(:,i) - LETbin_start(:,i));
    LETpoint(:,i) = (LETbin_start(:,i) + LETbin_end(:,i)) ./ 2;
    
    %     LETbin_start(:,i) = cell2mat(data(3:3 + n_bins-1, 1));
    %     LETbin_end(:,i) = cell2mat(data(3:3 + n_bins-1, 2));
    %     value(:,i) = cell2mat(data(3:3 + n_bins-1, 3));
    %
    %     LETbin_width(:,i) = abs(LETbin_end(:,i) - LETbin_start(:,i));
    %     LETpoint(:,i) = (LETbin_start(:,i) + LETbin_end(:,i)) ./ 2;
    
end

%% Estimate fluence and dose average LET spectra
for i = 1:n_detectors
    
    weightedsum = sum(LETbin_width(:,i) .* value(:,i));
    
    f(:,i) = value(:,i) ./ weightedsum;
    LET_f(i) = sum(LETpoint(:,i) .* f(:,i) .* LETbin_width(:,i));
    
    d(:,i) = LETpoint(:,i) .* f(:,i) ./ LET_f(i);
    LET_d(i) = sum(LETpoint(:,i) .* d(:,i) .* LETbin_width(:,i));
    
    sum(f(:,i) .* LETbin_width(:,i));
    sum(d(:,i) .* LETbin_width(:,i));
    
end

%% Plot LET spectrum at a selected depth

% str = {};
% counter = 0;
% 
% figure()
% hold on
% for i = [15 63 65 67 68] %  1:n_detectors % 27
%     
%     if ~any(isnan(f(:,i)))
%         counter = counter + 1;
%         plot(LETpoint(:,i), LETpoint(:,i) .* f(:,i), 'LineWidth', 1.0)
%         str{counter} = ...
%             ['Position ' num2str(counter) ' (' num2str(round(depth(1) + (i-1)*median(diff(depth)),4)) ' cm depth)'];
%     end
%     
% end
% xlabel('L (keV/\mum)')
% ylabel('L * f(L) (keV/\mum)')
% % legend(['depth = ' num2str(depth(1) + (i-1)*median(diff(depth))) ' cm'])
% legend(str)
% set(gca, 'FontSize', 14)

%% Interpolate
% depth_interp = depth(1):0.0005:depth(end);
% LET_d_interp = interp1(depth, LET_d, depth_interp, 'pchip');
% LET_f_interp = interp1(depth, LET_f, depth_interp, 'pchip');

%% Plot
% figure()
% plot(depth, LET_d, 'b-o', depth, LET_f, 'r-x') %, ...
%     %depth_interp, LET_d_interp, 'b-', depth_interp, LET_f_interp, 'r-')
% xlabel('Depth (cm)')
% ylabel('Mean LET (keV/\mum)')
% legend('LET_d', 'LET_f', 'Location', 'best')
%xlim([2.6 2.67])

% %% Write to XLSX file
% header = {'Depth (cm)', 'LET_d (keV/um)', 'LET_f (keV/um)'};
% xlswrite(destpath, header, 1, 'A1')
% xlswrite(destpath, depth.', 1, 'A2')
% xlswrite(destpath, LET_d.', 1, 'B2')
% xlswrite(destpath, LET_f.', 1, 'C2')

end