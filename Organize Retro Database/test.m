clear all;
close all;
clc;

%% Test
sourcepath = 'E:\Ny mappe'; % 'E:\UorganisertData'; % 'C:\Users\Delmon\Desktop\UnOrgDICOM'; % 'C:\Users\Delmon\Desktop\TESTUnOrgDICOM'; % 
destpath = 'E:\Organized data'; % 'C:\Users\Delmon\Desktop\OrgDICOM';
organizeDICOM(loadPatients(sourcepath), destpath);

%% Knalldata2
% Har overført T2- og warn.txt-filer til K for denne disken
% sourcepath = 'E:\Pat\New folder ';
% destpath = 'E:\Pat\T2files';

%% Knalldata4
%Har overført T2- og warn.txt-filer til K for denne disken
% sourcepath = 'E:\DUMP\New folder ';
% destpath = 'E:\T2 files'; % 'E:\DUMP\T2files';

%% Knalldata5
% Har overført T2- og warn.txt-filer til K for denne disken
% sourcepath = 'G:\Arpit\Pat\New folder ';
% destpath = 'G:\Arpit\Pat\T2files';

%%
% for i = 19
%     i
% 
%     % renameFoldersAndFiles([sourcepath '(' num2str(i) ')'], 8);
%     % organizeDICOM(loadPatients([path num2str(i)]));
%     
%     t2_tra_patients = getFiles(loadPatients( ...
%         [sourcepath '(' num2str(i) ')'] ), 'MR', 'T2', 'tra', 'ax');
%     copyFiles(t2_tra_patients, destpath);

% end

%% Edit DICOM files
% sourcepath = 'C:\Users\Delmon\Desktop\Test';
% getDICOMfolders(sourcepath)

%% Get all pretreatment T2 series
% input_xlsxpath = 'E:\Pasientliste_CervixRetrospektiv_IkkeSensitiv.xlsx';
% sourcepath = 'E:\TraT2pretreatFinal'; % 'E:\T2 files pretreat'; % 'E:\T2 files';
% destpath = 'E:\Temp'; % 'E:\T2 files pretreat';
% 
% data = readXLSXdocument(input_xlsxpath);
% CopyPreTreatPatients(sourcepath, destpath, data);
% fileID = fopen(fullfile(destpath, 'info2.txt'));
% C_header = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 1, ...
%    'Delimiter', ',');
% C_data = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', ...
%     'Delimiter', ',');
% C_data = reshape(vertcat(C_data{:}), length(C_data{1}), length(C_data));
% C_tot = [horzcat(C_header{:}); C_data];
% fclose(fileID);
% xlswrite(fullfile(destpath, 'info2.xlsx'), C_tot)
