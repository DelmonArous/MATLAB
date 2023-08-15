clear all
close all
fclose all
clc

sourcepath = 'C:\Users\Delmon\Desktop\Work\Cohort 1 plus Cohort 2';
%sourcepath = 'C:\Users\Delmon\Desktop\Work\TEST';

DWIROIpatients = readROIpatients(sourcepath);
DWIROIpatients = addADCinfo(sourcepath, DWIROIpatients);
DWIROIpatients = addfBVinfo(sourcepath, DWIROIpatients);

%% Reslice, resize T2w and T2* images and then generate R2* maps
% ResliceandInterpolate(sourcepath, 'T2', DWIROIpatients)
% ResliceandInterpolate(sourcepath, 'T2star', DWIROIpatients)

% generateT2starMap(sourcepath, 'T2starreslice_v2')
 
%%
R2starpatients = computeR2starROI(sourcepath, DWIROIpatients);
    
xlsxpath = 'C:\Users\Delmon\Desktop\Work\delmond.xlsx';
R2starpatients = R2starCorr(xlsxpath, R2starpatients);
  
% R2starPixelAnalysis(R2starpatients)
 
% T2wFatSignalPath = 'C:\Users\Delmon\Desktop\Work\T2 fat signal';
% R2starpatients = addT2wFatSI(T2wFatSignalPath, R2starpatients);
% R2starpatients = computeT2metric(sourcepath, R2starpatients);
%plotT2wmetricCorr(R2starpatients);

% checkAutoCorr(R2starpatients) 
% estimateHF_R2star(R2starpatients)
% R2starpatients = R2star_fBV_hypoxiaRegistration(R2starpatients);