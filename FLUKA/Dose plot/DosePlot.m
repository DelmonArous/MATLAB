clear all;
close all;
fclose('all');
clc;

path = 'C:\Users\Delmon\Desktop\16MeVProtonPhD';
sourcefilename = 'Dose.xlsx';

%% Read relevant data; differential fluence vs. energy
data = readXLSXdocument(fullfile(path,sourcefilename));

Dose = [];
bin_start = [];
bin_end = [];

Dose = cell2mat(data(2:end, 3));
bin_start = cell2mat(data(2:end, 1));
bin_end = cell2mat(data(2:end, 2));
depth = (bin_start(:) + bin_end(:)) ./ 2; % in cm
%depth = depth * 10; % in mm

maxVal = max(Dose(:)); 
BPind = find(Dose == maxVal); 
BPdepth = depth(BPind);

% Normalize to max
Dose_norm = (Dose ./ maxVal) * 100;
Dose_Gy = Dose .* 1.602176462*10^(-7); % in Gy

% Find R40, R60, R80, distal R80, distal R60 and distal R40
[~, indR40] = min(abs(Dose_norm(1:BPind) - 40));
R40 = depth(indR40);
[~, indR60] = min(abs(Dose_norm(1:BPind) - 60));
R60 = depth(indR60);
[~, indR80] = min(abs(Dose_norm(1:BPind) - 80));
R80 = depth(indR80);
[~, indDistalR80] = min(abs(Dose_norm(BPind:end) - 80));
DistalR80 = depth(BPind + indDistalR80 - 1);
[~, indDistalR60] = min(abs(Dose_norm(BPind:end) - 60));
DistalR60 = depth(BPind + indDistalR60 - 1);
[~, indDistalR40] = min(abs(Dose_norm(BPind:end) - 40));
DistalR40 = depth(BPind + indDistalR40 - 1);

%% Plot
figure();
hold on
plot(depth, Dose_norm)
plot(R40, Dose_norm(indR40), 'ro', R60, Dose_norm(indR60), 'ro', ...
    R80, Dose_norm(indR80), 'ro', BPdepth, Dose_norm(BPind), 'ro', ...
    DistalR80, Dose_norm(BPind + indDistalR80 - 1), 'ro', ...
    DistalR60, Dose_norm(BPind + indDistalR60 - 1), 'ro', ...
    DistalR40, Dose_norm(BPind + indDistalR40 - 1), 'ro')
text(R40, Dose_norm(indR40), 'R_{40}', 'HorizontalAlignment', 'right')
text(R60, Dose_norm(indR60), 'R_{60}', 'HorizontalAlignment', 'right')
text(R80, Dose_norm(indR80), 'R_{80}', 'HorizontalAlignment', 'right')
text(BPdepth, Dose_norm(BPind), 'BP MC', 'HorizontalAlignment', 'right')
text(DistalR80, Dose_norm(BPind + indDistalR80 - 1), 'R_{80,distal}', 'HorizontalAlignment', 'left')
text(DistalR60, Dose_norm(BPind + indDistalR60 - 1), 'R_{60,distal}', 'HorizontalAlignment', 'left')
text(DistalR40, Dose_norm(BPind + indDistalR40 - 1), 'R_{40,distal}', 'HorizontalAlignment', 'left')
hold off
legend(['16 MeV protons (BP ' num2str(round(BPdepth,2)) ' cm)'], ...
    'Location', 'best')
