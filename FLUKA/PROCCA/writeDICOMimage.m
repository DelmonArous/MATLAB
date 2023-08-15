function [varargout] = writeDICOMimage(img, path)

%% Initialize dicominfo structure with FileMetaInformationVersion
info.FileMetaInformationVersion = [0; 1];

%% Specify transfer syntax and implementation UIDs
info.TransferSyntaxUID          = '1.2.840.10008.1.2';
info.ImplementationClassUID     = '1.2.40.0.13.1.1';
info.ImplementationVersionName  = 'dcm4che-2.0';
info.SpecificCharacterSet       = 'ISO_IR 100';

% Specify class UID based on provided class
if isfield(img, 'MediaStorageSOPClassUID') && isfield(img, 'Modality')
    info.MediaStorageSOPClassUID    = img.MediaStorageSOPClassUID;
    info.Modality                   = img.Modality;
    
    % Otherwise, assume it is CT
else
    info.MediaStorageSOPClassUID    = '1.2.840.10008.5.1.4.1.1.2';
    info.Modality                   = 'CT';
end

% image.classUID = info.SOPClassUID;
% MediaStorageSOPClassUID
% image.studyUID = info.StudyInstanceUID;
% image.seriesUID = info.SeriesInstanceUID;
% image.frameRefUID = info.FrameOfReferenceUID;

%% Copy SOP class UID from Media
info.SOPClassUID = img.SOPClassUID; % info.MediaStorageSOPClassUID

%% If minimum image value is below zero, assume image is in HU and add 1024
if min(min(min(img.data))) < 0
    img.data = img.data + 1024;
end

%% Set intercept and slope based on modality
if strcmp(info.MediaStorageSOPClassUID, '1.2.840.10008.5.1.4.1.1.2') && ...
        isfield(img, 'RescaleIntercept') && isfield(img, 'RescaleSlope')
%     info.RescaleIntercept   = img.RescaleIntercept;
%     info.RescaleSlope       = img.RescaleSlope;
    info.RescaleIntercept   = -1024;
    info.RescaleSlope       = 1;
elseif strcmp(info.MediaStorageSOPClassUID, '1.2.840.10008.5.1.4.1.1.2') && ...
        ~isfield(img, 'RescaleIntercept') && ~isfield(img, 'RescaleSlope')
    info.RescaleIntercept   = -1024;
    info.RescaleSlope       = 1;
else
    info.RescaleIntercept   = 0;
    info.RescaleSlope       = 1;
end

%% Window settings for CT
if strcmp(info.MediaStorageSOPClassUID, '1.2.840.10008.5.1.4.1.1.2') && ...
        isfield(img, 'WindowCenter') && isfield(img, 'WindowWidth')
   info.WindowCenter = img.WindowCenter;
   info.WindowWidth  = img.WindowWidth;
end
if isfield(img, 'WindowCenterWidthExplanation')
    info.WindowCenterWidthExplanation = img.WindowCenterWidthExplanation;
end

%% Generate creation date/time
if isfield(img, 'timestamp')
    t = img.timestamp;
else
    t = now;
end

% Specify creation date/time
info.InstanceCreationDate = datestr(t, 'yyyymmdd');
info.InstanceCreationTime = datestr(t, 'HHMMSS');

% Specify acquisition date/time
info.AcquisitionDate = datestr(t, 'yyyymmdd');
info.AcquisitionTime = datestr(t, 'HHMMSS');

%% Specifty image type
info.ImageType = 'ORIGINAL/PRIMARY/AXIAL';

% Specify manufacturer, model, and software version
info.Manufacturer = ['MATLAB ', version];
info.ManufacturerModelName = 'WriteDICOMImage';
info.SoftwareVersion = '1.2';

%% Specify series description
if isfield(img, 'SeriesDescription')
    info.SeriesDescription = img.SeriesDescription;
else
    info.SeriesDescription = '';
end

%% Specify study description
if isfield(img, 'StudyDescription')
    info.StudyDescription = img.StudyDescription;
else
    info.StudyDescription = '';
end

%% Specify patient info (assume that if name isn't provided, nothing is provided)
if isfield(img, 'PatientName')
    
    % If name is a structure
    if isstruct(img.PatientName)
        
        % Add name from provided dicominfo
        info.PatientName = img.PatientName;
        
        % Otherwise, if name is a char array
    elseif ischar(img.PatientName)
        
        % Add name to family name
        info.PatientName.FamilyName = img.PatientName;
        
        % Otherwise, throw an error
    else
        error('Provided patient name is an unknown format');
    end
    
    % Specify patient ID
    if isfield(img, 'PatientID')
        info.PatientID = img.PatientID;
    else
        info.PatientID = '';
    end
    
    % Specify patient birthdate
    if isfield(img, 'PatientBirthDate')
        info.PatientBirthDate = img.PatientBirthDate;
    else
        info.PatientBirthDate = '';
    end
    
    % Specify patient sex
    if isfield(img, 'PatientSex')
        info.PatientSex = upper(img.PatientSex(1));
    else
        info.PatientSex = '';
    end
    
    % Specify patient age
    if isfield(img, 'PatientAge')
        info.PatientAge = img.PatientAge;
    else
        info.PatientAge = '';
    end
    
    % Otherwise, no data was provided and no UI exists to prompt user
else
    
    % Add generic data
    info.PatientName.FamilyName     = 'DOE';
    info.PatientName.GivenName      = 'J';
    info.PatientID                  = '00000000';
    info.PatientBirthDate           = '';
    info.PatientSex                 = '';
    info.PatientAge                 = '';
    
end

%% Specify slice thickness (in mm)
info.SliceThickness = img.width(3) * 10; % mm

%% Specify study UID
if isfield(img, 'StudyInstanceUID')
    info.StudyInstanceUID = img.StudyInstanceUID;
else
    info.StudyInstanceUID = dicomuid;
end

%% Specify series UID
if isfield(img, 'SeriesInstanceUID')
    info.SeriesInstanceUID = img.SeriesInstanceUID;
else
    info.SeriesInstanceUID = dicomuid;
end

%% Specify frame of reference UID
if isfield(img, 'FrameOfReferenceUID')
    info.FrameOfReferenceUID = img.FrameOfReferenceUID;
    % Otherwise, generate unique frame of reference ID
else
    info.FrameOfReferenceUID = dicomuid;
end

%% Specify patient position
if isfield(img, 'PatientPosition')
    info.PatientPosition = img.PatientPosition;
end

% Specify image orientation
if isfield(img, 'PatientPosition')
    
    % Set orientation based on patient position
    if strcmpi(img.PatientPosition, 'HFP')
        info.ImageOrientationPatient = [-1;0;0;0;-1;0];
    elseif strcmpi(img.PatientPosition, 'FFS')
        info.ImageOrientationPatient = [-1;0;0;0;1;0];
    elseif strcmpi(img.PatientPosition, 'FFP')
        info.ImageOrientationPatient = [1;0;0;0;-1;0];
    else
        info.ImageOrientationPatient = [1;0;0;0;1;0];
    end
    
    % Otherwise, assume standard (HFS) orientation
else
    info.ImageOrientationPatient = [1;0;0;0;1;0];
end

%% Specify IEC-X, IEC-Y and IEC-Z position, in mm
info.ImagePositionPatient(1) = img.start(1) * 10 * ...
    info.ImageOrientationPatient(1);
info.ImagePositionPatient(2) =  -(img.start(2) + img.width(2) * ...
    (size(img.data,2)-1)) * 10 * info.ImageOrientationPatient(5);
info.ImagePositionPatient(3) = -img.start(3) * 10;

%% Specify number of rows/columns and number of images
info.Rows = size(img.data, 1);
info.Columns = size(img.data, 2);
info.ImagesInAcquisition = size(img.data, 3);

%% Specify pixel spacing (in mm)
info.PixelSpacing = [img.width(1); img.width(2)] * 10;

%% Specify number of samples
info.SamplesPerPixel            = 1;
info.PhotometricInterpretation  = 'MONOCHROME2';

%% Specify bit information
info.BitsAllocated          = 16;
info.BitsStored             = 16;
info.HighBit                = 15;
info.PixelRepresentation    = 0;

%% Loop through images
for i = 1:size(img.data, 3)
    
    % Specify slice instance UID
    if isfield(img, 'SOPInstanceUID')
        info.MediaStorageSOPInstanceUID = img.SOPInstanceUID{i};
    else
        info.MediaStorageSOPInstanceUID = dicomuid;
    end
    
    % Specify return variable
    if nargout == 1
        varargout{1}{i} = info.MediaStorageSOPInstanceUID;
    end
    
    % Copy instance UID from media storage UID
    info.SOPInstanceUID = info.MediaStorageSOPInstanceUID;
    
    % Specify slice location in mm
    if isequal(info.ImageOrientationPatient, [1;0;0;0;1;0]) || ...
            isequal(info.ImageOrientationPatient, [-1;0;0;0;-1;0])
        % Slice position for head first
        info.SliceLocation = (img.start(3) + (i-1)*img.width(3)) * 10;
    else
        % Slice position for feet first
        info.SliceLocation = -(img.start(3) + (i-1)*img.width(3)) * 10;
    end
    
    % Update image position to slice location
    info.ImagePositionPatient(3) = -info.SliceLocation;
    
    % Add this slice's image instance number
    info.InstanceNumber = i;
    
    % Write DICOM file using dicomwrite()
    % uint16 for MR!!!
    status = dicomwrite(flip(rot90(uint16(img.data(:,:,i)), 3), 2), ...
        fullfile(path, sprintf('_%03i.dcm', i)), info, ...
        'CompressionMode', 'None', 'CreateMode', 'Copy', ...
        'Endian', 'ieee-le');
    
    % If status is not empty, break loop
    if ~isempty(status)
        break;
    end
    
end

end