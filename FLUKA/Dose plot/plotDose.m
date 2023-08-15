function [depth, Dose] = plotDose(sourcefile)

%% Read relevant data; differential fluence vs. energy
% data = readXLSXdocument(sourcefile);
A = importdata(sourcefile);

Dose = [];
bin_start = [];
bin_end = [];

% Dose = cell2mat(data(2:end, 3));
% bin_start = cell2mat(data(2:end, 1));
% bin_end = cell2mat(data(2:end, 2));
% depth = (bin_start(:) + bin_end(:)) ./ 2; % in cm
%depth = depth * 10; % in mm

Dose = A.data(1:end, 3);
bin_start = A.data(1:end, 1);
bin_end = A.data(1:end, 2);
depth = (bin_start(:) + bin_end(:)) ./ 2; % in cm

% maxVal = max(Dose(:));
% BPind = find(Dose == maxVal);

% Normalize to max
% Dose = Dose .* 1.602176462*10^(-7) * 10^9; % in nGy
% Dose = (Dose ./ maxVal) * 100;

% Find R40, R60, R80, distal R80, distal R60 and distal R40
% [~, indR40] = min(abs(Dose(1:BPind) - 40));
% [~, indR60] = min(abs(Dose(1:BPind) - 60));
% [~, indR80] = min(abs(Dose(1:BPind) - 80));
% [~, indDistalR80] = min(abs(Dose(BPind:end) - 80));
% [~, indDistalR60] = min(abs(Dose(BPind:end) - 60));
% [~, indDistalR40] = min(abs(Dose(BPind:end) - 40));
% R40 = depth(indR40);
% R60 = depth(indR60);
% R80 = depth(indR80);
% BPdepth = depth(BPind);
% DistalR80 = depth(BPind + indDistalR80 - 1);
% DistalR60 = depth(BPind + indDistalR60 - 1);
% DistalR40 = depth(BPind + indDistalR40 - 1);

%% Plot
% figure();
% hold on
% plot(depth, Dose, 'b')
% plot(R40, Dose(indR40), 'ro', R60, Dose(indR60), 'ro', ...
%     R80, Dose(indR80), 'ro', BPdepth, Dose(BPind), 'ro', ...
%     DistalR80, Dose(BPind + indDistalR80 - 1), 'ro', ...
%     DistalR60, Dose(BPind + indDistalR60 - 1), 'ro', ...
%     DistalR40, Dose(BPind + indDistalR40 - 1), 'ro')
% text(R40, Dose(indR40), 'R_{40}', 'HorizontalAlignment', 'right')
% text(R60, Dose(indR60), 'R_{60}', 'HorizontalAlignment', 'right')
% text(R80, Dose(indR80), 'R_{80}', 'HorizontalAlignment', 'right')
% text(BPdepth, Dose(BPind), 'BP MC', 'HorizontalAlignment', 'right')
% text(DistalR80, Dose(BPind + indDistalR80 - 1), 'R_{80,distal}', 'HorizontalAlignment', 'left')
% text(DistalR60, Dose(BPind + indDistalR60 - 1), 'R_{60,distal}', 'HorizontalAlignment', 'left')
% text(DistalR40, Dose(BPind + indDistalR40 - 1), 'R_{40,distal}', 'HorizontalAlignment', 'left')
% hold off
% xlabel('Depth (cm)')
% ylabel('Dose (%)')
% legend([str ' (BP ' num2str(round(BPdepth,2)) ' cm)'], 'Location', 'SouthWest')
% 
% figure();
% hold on
% plot(depth, Dose_Gy, 'b')
% plot(R40, Dose_Gy(indR40), 'ro', R60, Dose_Gy(indR60), 'ro', ...
%     R80, Dose_Gy(indR80), 'ro', BPdepth, Dose_Gy(BPind), 'ro', ...
%     DistalR80, Dose_Gy(BPind + indDistalR80 - 1), 'ro', ...
%     DistalR60, Dose_Gy(BPind + indDistalR60 - 1), 'ro', ...
%     DistalR40, Dose_Gy(BPind + indDistalR40 - 1), 'ro')
% text(R40, Dose_Gy(indR40), 'R_{40}', 'HorizontalAlignment', 'right')
% text(R60, Dose_Gy(indR60), 'R_{60}', 'HorizontalAlignment', 'right')
% text(R80, Dose_Gy(indR80), 'R_{80}', 'HorizontalAlignment', 'right')
% text(BPdepth, Dose_Gy(BPind), 'BP MC', 'HorizontalAlignment', 'right')
% text(DistalR80, Dose_Gy(BPind + indDistalR80 - 1), 'R_{80,distal}', ...
%     'HorizontalAlignment', 'left')
% text(DistalR60, Dose_Gy(BPind + indDistalR60 - 1), 'R_{60,distal}', ...
%     'HorizontalAlignment', 'left')
% text(DistalR40, Dose_Gy(BPind + indDistalR40 - 1), 'R_{40,distal}', ...
%     'HorizontalAlignment', 'left')
% hold off
% xlabel('Depth (cm)')
% ylabel('Dose (nGy/primary)')
% legend([str ' (BP ' num2str(round(BPdepth,2)) ' cm)'], 'Location', 'SouthWest')


end