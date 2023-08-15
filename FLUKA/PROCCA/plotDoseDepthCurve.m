clear all;
% close all;
fclose('all');
clc;

path = { ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu\DoseDepth_Xray_180kV_01mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.3 mm Cu\DoseDepth_Xray_180kV_03mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\200 kV 0.1 mm Cu\DoseDepth_Xray_200kV_01mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\200 kV 0.3 mm Cu\DoseDepth_Xray_200kV_03mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\220 kV 0.1 mm Cu\DoseDepth_Xray_220kV_01mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\220 kV 0.3 mm Cu\DoseDepth_Xray_220kV_03mmCu.xlsx'};

% path = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu - RCC\DoseDepth_Xray180kV_01mmCu.xlsx', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu - RPP\DoseDepth_Xray180kV_01mmCu.xlsx'};
 
str = {'180 kV, 0.1 mm Cu', '180 kV, 0.3 mm Cu', ...
    '200 kV, 0.1 mm Cu', '200 kV, 0.3 mm Cu', ...
    '220 kV, 0.1 mm Cu', '220 kV, 0.3 mm Cu'};
lineSpec = {'-r', '-b', '--r', '--b', '-.r', '-.b'};

% str = {'180 kV, 0.1 mm Cu, Cylinder', '180 kV, 0.1 mm Cu, Box'};
% lineSpec = {'-r', '--r'};

figure();
hold on
for i = 1:length(path)
   
    [depth, dose] = readDoseDepth(path{i});
    plot(depth, dose, lineSpec{i})
    
end
hold off
ylim([-inf inf])
xlabel('Depth (cm)')
ylabel('Dose (nGy per primary)')
legend(str, 'Location', 'NorthEast')
set(gca, 'FontSize', 14)
