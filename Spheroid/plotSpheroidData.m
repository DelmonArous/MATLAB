clear all;
close all
clc

%% Read T47D spheroid data
sourcepath = 'C:\Users\Delmon\Dropbox\Jobb\Matlab\Spheroid\T47D\T47D\Excel data';
fileList = getAllFiles(sourcepath);

spheroid_name = [];
outer_area_mm = [];
inner_area_mm = [];
growth = [7 9 11 13 17 24];

for i = 1:length(fileList)
    
    [p, n, e] = fileparts(fileList{i});
    data = readXLSXdocument(fullfile(p, [n e]));
   
    spheroid_name = [spheroid_name data(2:end,1)];    
    outer_area_mm = [outer_area_mm cell2mat(data(2:end, 2))];
    inner_area_mm = [inner_area_mm cell2mat(data(2:end, 5))];
    
end

%% 
a = [1 11 20 30 40 50];
b = [10 19 29 39 49 59];

tempWell_mean_outer_area_mm = [];
tempWell_mean_inner_area_mm = [];
n_wells = [];

for i = 1:length(fileList)
    
    for j = 1:6
        
        n_wells(j) = length(a(j):b(j));
        tempWell_outer_area_mm = rmmissing(outer_area_mm(a(j):b(j), i));
        tempWell_inner_area_mm = rmmissing(inner_area_mm(a(j):b(j), i));
  
        % Mean and Standard Deviation
        tempWell_mean_outer_area_mm(j,i) = mean(tempWell_outer_area_mm);
        tempWell_mean_inner_area_mm(j,i) = mean(tempWell_inner_area_mm);
        SD_outer(j,i) = std(tempWell_outer_area_mm);
        SD_inner(j,i) = std(tempWell_inner_area_mm);

    end

end

%% Plot
linespec = {'-+', '--o', ':p', '-.d', '-x', '--s'};

figure()
hold on
for i = 1:6
   
   yyaxis left
   errorbar(growth, tempWell_mean_outer_area_mm(i,:), SD_outer(i,:), linespec{i})
   ylabel('Spheroid Cross-Sectional Area (mm^2)')
   yyaxis right
   errorbar(growth, tempWell_mean_inner_area_mm(i,:), SD_inner(i,:), linespec{i})
   ylabel('Necrotic Center Cross-Sectional Area (mm^2)')
 
end
lgd = legend(['Well B, n=' num2str(n_wells(1))], ...
    ['Well C, n=' num2str(n_wells(2))], ...
    ['Well D, n=' num2str(n_wells(3))], ...
    ['Well E, n=' num2str(n_wells(4))], ...
    ['Well F, n=' num2str(n_wells(5))], ...
    ['Well G, n=' num2str(n_wells(6))], 'Location', 'NorthWest');
title(lgd, ['n=' num2str(size(spheroid_name,1)), ' T47D Spheroids'])
xlim([growth(1)-1 growth(end)+1])
xticks(growth)
xlabel('Growth (days)')
