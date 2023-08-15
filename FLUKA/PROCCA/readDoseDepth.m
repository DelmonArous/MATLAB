function [depth, Dose] = readDoseDepth(path)

%% Read relevant data

bin_start = []; bin_end = [];
Dose = [];

data = readXLSXdocument(path);

bin_start = cell2mat(data(2:end, 1));
bin_end = cell2mat(data(2:end, 2));
Dose = cell2mat(data(2:end, 3)) .* 1.602176462*10^(-7) .* 10^9; % in nGy

depth = (bin_start + bin_end) ./ 2 ; 
% depth = depth .* 10; % in mm

end