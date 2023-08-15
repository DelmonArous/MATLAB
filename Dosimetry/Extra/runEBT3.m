clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path = 'C:\Users\delmo\Dropbox\Jobb\Matlab\Dosimetry\EBT3 Olga\New Folder\Data';

%% Variables

% Scanning resolution (in dpi), bit depth and ROI size (in mm)
dpi = 300;
bit = 48;
ROI_size_mm = [4 4];

% EBT3 film dose (in Gy)
dose = [10.0 10.0 10.0 10.0 5.0 5.0 5.0 5.0 7.0 7.0 7.0 7.0];

% Define EBT3 film x- and y-coordinate range (in pixels)
x_range_px = [1  580];
y_range_px = [35 613];

%% Get EBT3 variable structure
EBT3struct = getEBT3struct(path, dpi, bit, ROI_size_mm, ...
    x_range_px, y_range_px);




