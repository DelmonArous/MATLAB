clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path = 'C:\Users\Delmon Arous\Desktop\EBT3';

%% Variables

% Scanning resolution (in dpi), bit depth and ROI size (in mm)
dpi = 300; % [1200 300 1200 300];
bit = 48;
ROI_size_mm = [4 4];

% EBT3 film dose (in Gy)
% dose = {[10.0 10.0 10.0 10.0 5.0 5.0 5.0 5.0 7.0 7.0 7.0 7.0], ...
%     [10.0 10.0 10.0 10.0 5.0 5.0 5.0 5.0 7.0 7.0 7.0 7.0], ...
%     [5.0 5.0 5.0 5.0], [5.0 5.0 5.0 5.0]};
doses = {[  ]}; 

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = {[1 2280], [1  570], [1 2000], [1 500]};
y_range_px = {[100 2340], [20 580], [50 2850], [10 710]};

% Plot spesifications
channel     = {'red', 'green', 'blue', 'gray'};
markerspec  = {'ro', 'gx', 'bs', 'kd'};
linespec    = {'r-', 'g-', 'b-', 'k-'};
lgdstr      = {'Red channel', 'Green channel', 'Blue channel', ...
    'Grayscale image'};

%%

EBT3films = getEBT3struct(path, dpi, bit, ...
    ROI_size_mm, x_range_px, y_range_px);



%%

% folderlist = getAllFolders(path);

% for i = 3:4 % :length(folderlist)
    
%     if i == 3
%         EBT3films_1200dpi = getEBT3struct(folderlist{i}, dpi(i), bit, ...
%         ROI_size_mm, x_range_px{i}, y_range_px{i});
%     else
%         EBT3films_300dpi = getEBT3struct(folderlist{i}, dpi(i), bit, ...
%         ROI_size_mm, x_range_px{i}, y_range_px{i});
%     end
    
%     EBT3films = getEBT3struct(folderlist{i}, dpi(i), bit, ...
%         ROI_size_mm, x_range_px{i}, y_range_px{i});
    
%     figure();
%     hold on
%     for j = 1:length(channel)
%         
%         PVvec       = [];
%         ODvec       = [];
%         
%         for k = 1:length(EBT3films)
%             
%             PVvec = [PVvec EBT3films{k}.(channel{j}).PV];
%             ODvec = [ODvec EBT3films{k}.(channel{j}).OD];
%             
%         end
%         
%         h(j) = plot(dose{i}, PVvec, markerspec{j});
%         
%     end
%     
%     xlabel('Dose (Gy)')
%     ylabel('\it{PV}')
%     xlim([min(dose{i})-0.5 max(dose{i})+0.5])
%     % ylim([11 12])
%     legend([h(1), h(2), h(3), h(4)], ...
%         'Data, red channel', 'Data, green channel', 'Data, blue channel', ...
%         'Data, grayscale image', 'Location', 'NorthEast')
%     grid on
%     title(['EBT3 scan ' num2str(dpi(i)) ' dpi'])
%     set(gca, 'FontSize', 16)
    
% end

%%

for k = 1:length(EBT3films_1200dpi) 
            
%     x_1200dpi = [size(EBT3films_1200dpi{k}.rgb.img,2)/2 ...
%         size(EBT3films_1200dpi{k}.rgb.img,2)/2];
%     y_1200dpi = [0 size(EBT3films_1200dpi{k}.rgb.img,1)];
%     c_1200dpi = improfile(EBT3films_1200dpi{k}.rgb.img, ...
%         x_1200dpi, y_1200dpi);
%     
%     x_300dpi = [size(EBT3films_300dpi{k}.rgb.img,2)/2 ...
%         size(EBT3films_300dpi{k}.rgb.img,2)/2];
%     y_300dpi = [0 size(EBT3films_300dpi{k}.rgb.img,1)];
%     c_300dpi = improfile(EBT3films_300dpi{k}.rgb.img, ...
%         x_300dpi, y_300dpi);
   
    y_temp = linspace(EBT3films_1200dpi{i}.y(1), ...
        EBT3films_1200dpi{i}.y(2), length(EBT3films_300dpi{i}.c(:,1,1)))';

    figure();
    hold on
    subplot(2,1,1)
    imshow(EBT3films_1200dpi{k}.rgb.img)
    hold on
    plot(EBT3films_1200dpi{k}.x, EBT3films_1200dpi{k}.y, 'r')
    subplot(2,1,2)
    plot(EBT3films_1200dpi{k}.c(:,1,1), '-r')
    plot(y_temp, EBT3films_300dpi{k}.c(:,1,1), '--r')
    hold on
    plot(EBT3films_1200dpi{k}.c(:,1,2), '-g')
    plot(y_temp, EBT3films_300dpi{k}.c(:,1,2), '--g')
    plot(EBT3films_1200dpi{k}.c(:,1,3), '-b')
    plot(y_temp, EBT3films_300dpi{k}.c(:,1,3), '--b')
    hold off
%     title(strrep(EBT3films_1200dpi{k}.filename, '_', '-'))
    xlabel('y')
    ylabel('\it{PV}')
    set(gca, 'FontSize', 14)

end
