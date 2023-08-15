function [LET_f, LET_d] = averageLET(sourcepath)

filelist    = getAllFiles(sourcepath);
n_detectors = length(filelist); % number of detectors (LET measurements)

%% Read data

% Initialization
LETbin_start = []; LETbin_end = []; LETbin_width = []; LETpoint = [];
value = [];
f = []; d = []; LET_f = []; LET_d = [];

formatSpec = '%f %f %f %f';
dataarray_size = [4 Inf];

% Import data into arrays
for i = 1:n_detectors

    %data = readXLSXdocument(fullfile(sourcepath, ['LET' num2str(i-1) '.xlsx']));
    %     fileID = fopen(fullfile(sourcepath, ['PhotonGrid_16MV_' num2str(i+29) '_tab.lis']), 'r');
    fileID = fopen(filelist{i}, 'r');
    fgetl(fileID); fgetl(fileID); % skip the first two lines
    dataarray = fscanf(fileID, formatSpec, dataarray_size);
    dataarray =  dataarray';

    LETbin_start(:,i)   = dataarray(:,1);
    LETbin_end(:,i)     = dataarray(:,2);
    value(:,i)          = dataarray(:,3);

    LETbin_width(:,i)   = abs(LETbin_end(:,i) - LETbin_start(:,i));
    LETpoint(:,i)       = (LETbin_start(:,i) + LETbin_end(:,i)) ./ 2;

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

    % Fluence-weighted LET
    f(:,i)      = value(:,i) ./ weightedsum;
    LET_f(i)    = sum(LETpoint(:,i) .* f(:,i) .* LETbin_width(:,i));

    % Dose-weighted LET
    d(:,i)      = LETpoint(:,i) .* f(:,i) ./ LET_f(i);
    LET_d(i)    = sum(LETpoint(:,i) .* d(:,i) .* LETbin_width(:,i));

    sum(f(:,i) .* LETbin_width(:,i))
    sum(d(:,i) .* LETbin_width(:,i))

end

end