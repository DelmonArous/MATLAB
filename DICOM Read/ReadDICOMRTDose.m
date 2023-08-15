function [dose] = ReadDICOMRTDose(varargin)

% The function loads a DICOM RTDose object into a MATLAB structure that
% can be used for manipulation with other functions in this library.

% The following variables are required for proper execution: 
%   varargin{1}: string containing the path to the DICOM files
%   varargin{2}: string containing the DICOM RTDose file to be loaded

% The following variables are returned upon succesful completion:
%   dose: structure containing the image data, dimensions, width, start 
%       coordinates, and key DICOM header values. The data is a three 
%       dimensional array of dose values in the units specified in the 
%       DICOM header, while the dimensions, width, and start fields are 
%       three element vectors. The DICOM header values are returned as 
%       strings.

%% Read DICOM header info from file
info = dicominfo(fullfile(varargin{1}, varargin{2}));

%% Read Unique Identifier (UID) and patient characteristics
dose.SOPClassUID = info.SOPClassUID;
dose.StudyInstanceUID = info.StudyInstanceUID;
dose.SeriesInstanceUID = info.SeriesInstanceUID;
dose.FrameOfReferenceUID = info.FrameOfReferenceUID;
dose.ImagePositionPatient = info.ImagePositionPatient; % start voxel coordinate

if isfield(info, 'PatientName')
    dose.PatientName = info.PatientName;
end
if isfield(info, 'PatientID')
    dose.PatientID = info.PatientID;
end
if isfield(info, 'PatientBirthDate')
    dose.PatientBirthDate = info.PatientBirthDate;
end
if isfield(info, 'PatientSex')
    dose.PatientSex = info.PatientSex;
end
if isfield(info, 'PatientAge')
    dose.PatientAge = info.PatientAge;
end

%% If patient orientation does not exist, assume HFS
%if ~isfield(info, 'ImageOrientationPatient') 
    info.ImageOrientationPatient = [1;0;0;0;1;0];
%end

%% Check if slice positions vary more than 1% 
widthsdiff = diff(info.GridFrameOffsetVector); % grid offset differences in mm
if ( abs(max(widthsdiff) - min(widthsdiff))/mean(widthsdiff) > 0.01 )
    error(['Error occured. \n', ...
        'Slice position vary more than 1% because of degree of unequal grid spacing.']);
end

% Store pixel and voxel dimensions in mm
dose.width = info.PixelSpacing; % X and Z voxel dimensions in cm
dose.width(3) = abs(mean(widthsdiff)); % mean IEC-Y width in cm

%% If GridFrameOffsetVector is descending, read in dose data and convert to double
S = info.DoseGridScaling; % D_i = V_i*S
if (info.GridFrameOffsetVector(end) < info.GridFrameOffsetVector(1))
    dose.values = flip(rot90(double(squeeze(dicomread(info)))))*S;
else % Else GridFrameOffsetVector is ascending
    dose.values = flip(flip(rot90(double(squeeze(dicomread(info))))), 3)*S;
end

% Resolution of the dose image
dose.dimensions = size(dose.values);

%% Start IEC-X voxel coordinate in mm
dose.start(1) = info.ImagePositionPatient(1) * ...
    info.ImageOrientationPatient(1);

%% Start IEC-Z voxel coordinate (inverted value) in mm 
dose.start(2) = -(info.ImagePositionPatient(2) * ...
    info.ImageOrientationPatient(5) + dose.width(2) * ...
    (dose.dimensions(2) - 1));

%% Start IEC-Y voxel coordinate in cm based on patient position
if isequal(info.ImageOrientationPatient, [1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [-1;0;0;0;-1;0]) 
    
    % Store patient position, head first (HF)
    if (info.ImageOrientationPatient(5) == 1)
        dose.position = 'HFS';
    elseif (info.ImageOrientationPatient(5) == -1)
        dose.position = 'HFP';
    end
    
    dose.start(3) = -max(info.ImagePositionPatient(3) + ...
        info.GridFrameOffsetVector);
    
else
    
    % Store patient position, feet first (FF)
    if (info.ImageOrientationPatient(5) == 1)
        dose.position = 'FFS';
    elseif (info.ImageOrientationPatient(5) == -1)
        dose.position = 'FFP';
    end
    
    dose.start(3) = max(info.ImagePositionPatient(3) + ...
        info.GridFrameOffsetVector);
    
end

%% Clear temporary variables
clear info widthsdiff;

end