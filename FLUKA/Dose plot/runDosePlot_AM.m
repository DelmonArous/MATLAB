clear all;
close all;
fclose('all');
clc;

% sourcepath = {'C:\Users\delmo\Desktop\AM paper\Dose\Proton 14.5 MeV\XLSX', ...
%     'C:\Users\delmo\Desktop\AM paper\Dose\Proton 15.5 MeV\XLSX', ...
%     'C:\Users\delmo\Desktop\AM paper\Dose\Proton 16.5 MeV\XLSX'};
% C:\Users\delmo\Desktop\AM paper\Fluence\XLSX
% sourcepath = {'C:\Users\delmo\Desktop\AM paper\Dose\Lid Parafilm\XLSX'};
% sourcepath = {'C:\Users\delmo\Desktop\AM paper\AM data\AMproton_V2\XLSX v4'};
% sourcepath = {'C:\Users\delmo\Desktop\AM paper\AM data\AMproton_V1\XLSX\Lid', ...
%     'C:\Users\delmo\Desktop\AM paper\AM data\AMproton_V1\XLSX\Parafilm'};
% sourcepath = {'C:\Users\delmo\Desktop\AM paper\Dose\Proton 15.2 MeV Water target\XLSX'};
% sourcepath = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\AM paper\Dose\Proton 15.22 MeV 103.1 cm Dist\Lid', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\AM paper\Dose\Proton 15.22 MeV 103.1 cm Dist\Parafilm'};
sourcepath = {'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\AM paper\Dose\Proton 15.22 MeV\v2'};

% str_E = {'14.5 MeV protons', '15.5 MeV protons', '16.5 MeV protons'};
% str_dist = {'30 cm', '40 cm', '50 cm', '60 cm', '70 cm', '80 cm', ...
%     '90 cm', '100 cm', '110 cm', '120 cm', '130 cm'};
% str = {'15.0 MeV - Lid', '15.0 MeV - Parafilm', ...
%     '15.5 MeV - Lid', '15.5 MeV - Parafilm', ...
%     '16.0 MeV - Lid', '16.0 MeV - Parafilm', ...
%     '16.5 MeV - Lid', '16.5 MeV - Parafilm'};
% str = {'70 cm - Lid', '80 cm - Lid', '90 cm - Lid', ...
%     '70 cm - Parafilm', '80 cm - Parafilm', '90 cm - Parafilm'};
% str = {'0.100 mm', '0.980 mm', '0.985 mm', '0.990 mm', '0.995 mm', ...
%     '1.000 mm', '1.005 mm', '1.025 mm', '1.055 mm', '1.080 mm'};
% str = {'15.0 MeV', '15.1 MeV', '15.2 MeV', '15.3 MeV', '15.4 MeV', ...
%     '15.5 MeV'};
% lgdtitle = {'Lid, 81.3 cm', 'Parafilm, 81.3 cm'};
% str = {'15.10 MeV', '15.11 MeV', '15.12 MeV', '15.13 MeV', '15.14 MeV', ...
%     '15.15 MeV', '15.16 MeV', '15.17 MeV', '15.18 MeV', '15.19 MeV', ...
%     '15.20 MeV', '15.21 MeV', '15.22 MeV', '15.23 MeV', '15.24 MeV', ...
%     '15.25 MeV', '15.26 MeV', '15.27 MeV', '15.28 MeV', '15.29 MeV', ...
%     '15.30 MeV'};
str = {'Pos 1: 103.1 cm (Parafilm)', 'Pos 2: 95.2 cm (Lid)', ...
    'Pos 3: 98.5 cm (Lid)', 'Pos 4: 101.1 cm (Lid)', 'Pos 5: 103.1 cm (Lid)'};
% lgdtitle = {'Lid, pos. 103.1 cm', 'Parafilm, pos. 103.1 cm'};
lgdtitle = {'IC position'};

lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y', ...
    '--r', '--g', '--c', '--m', '--k', '--y', ...
    ':r', ':g', ':c', ':m', ':k', ':y', ...
    '-.r', '-.g', '-.c', '-.m', '-.k', '-.y'};

% corr_dist = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
% corr_dist = [54.735 55.687 54.735 55.687 54.735 55.687 54.735 55.687];
% corr_dist = [70-0.755 80-0.755 90-0.755 70-0.013 80-0.013 90-0.013];
% corr_dist = 88.0-1.445;
% corr_dist = [81.3-0.965 81.3-0.013]; % [Lid Parafilm]
corr_dist = [81.3-0.013 73.4-0.965 76.7-0.965 79.3-0.965 81.3-0.965]; % [Parafilm Lid Lid Lid Lid]

%%
for i = 1:length(sourcepath)
    
    fileList = getAllFiles(sourcepath{i});
    
    figure();
    hold on
    for j = 1:length(fileList)
                
        [filepath, filename, ext] = fileparts(fileList{j});
        [depth, Dose] = plotDose(fullfile(filepath, [filename ext]));

        % Convert to Gy
%         Dose = Dose .* 1.602176462*10^(-7) * 10^9; % in nGy
        % Normalize to max
        Dose = (Dose ./ max(Dose(:))) * 100;
        % Correct file depth
        depth = depth - corr_dist(j); % corr_dist(i);
        
        % Plot
        plot(depth, Dose, lineSpec{j}, 'LineWidth', 1)

    end
    hold off
    xlim([0 inf])
    ylim([0 inf])
    xlabel('Depth (cm)')
    ylabel('Dose (%)') % nGy per primary, particles/cm^2 per primary
    lgd = legend(str, 'Location', 'NorthEast');
    title(lgd, lgdtitle{i})
%     lgd = legend(str_dist, 'Location', 'NorthEast');
%     title(lgd, str_E{i})
    set(gca, 'FontSize', 14)
    
end
