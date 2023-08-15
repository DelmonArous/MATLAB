function [] = plotFluence(sourcepath, str)

%% Get all files in directory

fileList = getAllFiles(sourcepath);

%% Read relevant data; differential fluence vs. energy

Fluence = [];
bin_start = [];
bin_end = [];

for i = 1:length(fileList)
    
    fileList{i}
    
    data = readXLSXdocument(fileList{i});
    bin_start(:,i) = cell2mat(data(2:end, 1));
    bin_end(:,i) = cell2mat(data(2:end, 2));
    Fluence(:,i) = cell2mat(data(2:end, 3));
    depth(:,i) = (bin_start(:,i) + bin_end(:,i)) ./ 2;
    %depth(:,i) = depth(:,i) .* 10; % in mm
    
end

%% Plot
lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y'};
figure();
hold on
for i = 1:length(fileList)
    plot(depth(:,i), Fluence(:,i), lineSpec{i})
end
hold off
ylim([0 inf])
xlabel('Depth (cm)')
ylabel('Fluence (particles/cm^2 per primary)')
legend(str, 'Location', 'NorthEast')

for i = 1:length(fileList)
    figure();
    plot(depth(:,i), Fluence(:,i), lineSpec{i})
    xlabel('Depth (cm)')
    ylabel('Fluence (particles/cm^2 per primary)')
    legend(str{i}, 'Location', 'NorthEast')
end

end