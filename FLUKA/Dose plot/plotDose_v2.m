function [] = plotDose_v2(sourcepath, str)

%% Get all files in directory

fileList = getAllFiles(sourcepath);

%% Read relevant data

Dose = [];
bin_start = [];
bin_end = [];

for i = 1:length(fileList)
    
    fileList{i}
    
    data = readXLSXdocument(fileList{i});
    bin_start(:,i) = cell2mat(data(2:end, 1));
    bin_end(:,i) = cell2mat(data(2:end, 2));
    Dose(:,i) = cell2mat(data(2:end, 3)) .* 1.602176462*10^(-7) .* 10^9; % in nGy
    depth(:,i) = (bin_start(:,i) + bin_end(:,i)) ./ 2;
    %depth(:,i) = depth(:,i) .* 10; % in mm
    
end

%% Plot
lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y'};
figure();
hold on
for i = 1:length(fileList)
    plot(depth(:,i), Dose(:,i), lineSpec{i})
end
hold off
ylim([0 inf])
xlabel('Depth (cm)')
ylabel('Dose (nGy per primary)')
legend(str, 'Location', 'NorthEast')

for i = 1:length(fileList)
    figure();
    plot(depth(:,i), Dose(:,i), lineSpec{i})
    xlabel('Depth (cm)')
    ylabel('Dose (nGy per primary)')
    legend(str{i}, 'Location', 'NorthEast')
end

end