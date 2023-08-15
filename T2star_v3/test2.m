clear all
close all
clc

sourcepath = 'C:\Users\Delmon\Desktop\Work\ReslicedData';

DWIROIpatients = readROIpatients(sourcepath);
DWIROIpatients = addADCinfo(sourcepath, DWIROIpatients);
DWIROIpatients = addfBVinfo(sourcepath, DWIROIpatients);
  
%% Run for good interp patients - reslice and then generate R2* maps
% ResliceandInterpolate(sourcepath, 'T2', DWIROIpatients)
% ResliceandInterpolate(sourcepath, 'T2star', DWIROIpatients)
% 
% generateT2starMap(sourcepath, 'T2starreslice_v2')

%% Run only for bad interp patients - generate R2* maps and then reslice
%%% Put all such patients in a separate folder before execution!
% sourcepath = '.......';
% 
% generateT2starMap(sourcepath, 'T2star')
%
% ResliceandInterpolate(sourcepath, 'T2', DWIROIpatients)
% ResliceandInterpolate(sourcepath, 'R2star', DWIROIpatients)

fclose all;
fclose('all');