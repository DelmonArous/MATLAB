clear all;
close all;
clc;

sourcepath = 'C:\Users\Delmon\Dropbox\Jobb\Matlab\Signal denoising\RB.txt';

fileID = fopen(fullfile(sourcepath));
C = textscan(fileID, '%s %s', 'Delimiter', ',');
fclose(fileID);

%C = reshape(vertcat(C{:}), length(C{1}), length(C));

x = cell2mat(C{1});
y = cell2mat(C{2}.');


%data = readXLSXdocument(input_xlsxpath);
%C_data = reshape(vertcat(C_data{:}), length(C_data{1}), length(C_data));
%C_tot = [horzcat(C_header{:}); C_data];
