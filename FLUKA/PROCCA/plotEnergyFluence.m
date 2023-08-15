clear all;
% close all;
fclose('all');
clc;

path = { ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.1 mm Cu\EnergyFluence_Xray_180kV_01mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\180 kV 0.3 mm Cu\EnergyFluence_Xray_180kV_03mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\200 kV 0.1 mm Cu\EnergyFluence_Xray_200kV_01mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\200 kV 0.3 mm Cu\EnergyFluence_Xray_200kV_03mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\220 kV 0.1 mm Cu\EnergyFluence_Xray_220kV_01mmCu.xlsx', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Cylinder Data\220 kV 0.3 mm Cu\EnergyFluence_Xray_220kV_03mmCu.xlsx'};

str = {'180 kV, 0.1 mm Cu', '180 kV, 0.3 mm Cu', ...
    '200 kV, 0.1 mm Cu', '200 kV, 0.3 mm Cu', ...
    '220 kV, 0.1 mm Cu', '220 kV, 0.3 mm Cu'};
lineSpec = {'-r', '-b', '--r', '--b', '-.r', '-.b'};

ndetectors = 1;
nbins = 400;

figure();
hold on
for i = 1:length(path)
   
    [E, dPhi] = readEnergyFluence(path{i}, ndetectors, nbins);
    plot(E(:,1), dPhi(:,1), lineSpec{i})
    
end
hold off
xlim([0 220])
ylim([0 inf])
xlabel('Energy (kV)')
ylabel('Relative fluence (particles/cm^2/kV)') % Relative fluence (1/cm^2/kV)
legend(str, 'Location', 'NorthEast')
set(gca, 'FontSize', 14)
