function img = readDICOMImages(path, filenames)

% The following variables are required for proper execution: 
%   path: string containing the path to the DICOM files
%   names: cell array of strings containing all files to be loaded
%
% The following variables are returned upon succesful completion:
%   image: structure containing the image data, dimensions, width, type,
%       start coordinates, and key DICOM header values. The data is a three 
%       dimensional array of image values, while the dimensions, width, and 
%       start fields are three element vectors.  The DICOM header values 
%       are returned as a strings.
%
% Below is an example of how this function is used:
%
%   path = '/path/to/files/';
%   names = {
%       '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.1.dcm'
%       '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.2.dcm'
%       '2.16.840.1.114362.1.5.1.0.101218.5981035325.299641582.274.3.dcm'
%   };
%   image = LoadDICOMImages(path, names);

%% Check if MATLAB can find dicominfo (Image Processing Toolbox)
if exist('dicominfo', 'file') ~= 2
    
    % If not, throw an error
    error(['The Image Processing Toolbox cannot be found and is ', ...
        'required by this function.']);
    
end

% Execute in try/catch statement
try

%% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    
    % Start waitbar
    progress = waitbar(0, 'Loading DICOM images');
    
end

%% Initialize empty variables for the UIDs, patient demographics, and image dimensions
img.SOPClassUID           = '';
img.StudyInstanceUID      = '';
img.SeriesInstanceUID     = '';
img.FrameOfReferenceUID   = '';
img.SOPInstanceUID        = cell(0);
img.PatientName           = '';
img.PatientID             = '';
img.PatientBirthDate      = '';
img.PatientSex            = '';
img.PatientAge            = '';

%% Loop through each file in filename list

% Initialize empty 3D array for images and vector of slice locations
% (the data may not be loaded in correct order; these will be used to
% re-sort the slices later)
images = [];
sliceLocations = [];

for i = 1:length(filenames)
    
    % Update waitbar
    if exist('progress', 'var') && ishandle(progress)
        waitbar(i/(length(filenames)+2), progress, ...
            sprintf('Loading DICOM images (%i/%i)', i, length(filenames)));
    end
    
    % Attempt to load each file using dicominfo
    try
        
        % If dicominfo is successful, store the header information
        info = dicominfo(fullfile(path, filenames{i}));
        
    catch
        
        % Otherwise, the file is either corrupt or not a real DICOM
        % file, so warn user
        warning(['File ', filenames{i}, ' is not a valid DICOM image ', ...
            'and was skipped']);

        % Then, automatically skip to next file in directory 
        continue
    end 
    
    % If this is the first DICOM image (and the class UID
    % have not yet been set
    if strcmp(img.SOPClassUID, '') % classUID

        % Store the UIDs, patient demographics
        img.SOPClassUID           = info.SOPClassUID;
        img.StudyInstanceUID      = info.StudyInstanceUID;
        img.SeriesInstanceUID     = info.SeriesInstanceUID;
        img.FrameOfReferenceUID   = info.FrameOfReferenceUID;
        
        if isfield(info, 'PatientName')
            img.PatientName = info.PatientName;
        end
        if isfield(info, 'PatientID')
            img.PatientID = info.PatientID;
        end
        if isfield(info, 'PatientBirthDate')
            img.PatientBirthDate = info.PatientBirthDate;
        end
        if isfield(info, 'PatientSex')
            img.PatientSex = info.PatientSex;
        end
        if isfield(info, 'PatientAge')
            img.PatientAge = info.PatientAge;
        end
        
    % Otherwise, if this file's study UID does not match the others,
    % multiple DICOM studies may be present in the same folder (not
    % currently supported)
    elseif ~strcmp(img.StudyInstanceUID, info.StudyInstanceUID)
        error(['Multiple DICOM Study Instance UIDs were found in ', ...
            'this list.  Please select only one study.']);
        
    % Otherwise, if this file's series UID does not match the others,
    % multiple DICOM series may be present in the same folder (not
    % currently supported)
    elseif ~strcmp(img.SeriesInstanceUID, info.SeriesInstanceUID) 
        
        error(['Multiple DICOM Series Instance UIDs were found in ', ...
            'this list.  Please select only one series.']);
        
    end
    
    % Append this slice's instance UID
    img.SOPInstanceUID{length(img.SOPInstanceUID) + 1} = ...
        info.SOPInstanceUID;
    
    % Append this slice's location to the sliceLocations vector
    sliceLocations(length(sliceLocations) + 1) = ...
        info.ImagePositionPatient(3);
    
    % Append this slice's image data to the images array
    images(size(images,1)+1,:,:) = dicomread(info);
        
end

%% Update waitbar
if exist('progress', 'var') && ishandle(progress)
    waitbar((length(filenames)+1)/(length(filenames)+2), progress, ...
        'Processing images');
end

%% Verify info is set (if not no valid DICOM files were found)
if ~exist('info', 'var')
    
    error(['No valid DICOM images were found in the provided directory.', ...
        'Try running dicominfo() on the target file to see what error ', ...
        'occurred.']);
    
end

%% Set image type based on DICOM header modality tag
if isfield(info, 'Modality')
    img.Modality = info.Modality;
end

%% Retrieve start voxel IEC-X coordinate from DICOM header (in cm)
img.start(1) = info.ImagePositionPatient(1) / 10 * ...
    info.ImageOrientationPatient(1);

%% Adjust IEC-Z to inverted value (in cm)
img.start(2) = -(info.ImagePositionPatient(2) * ...
    info.ImageOrientationPatient(5) + info.PixelSpacing(2) * ...
    (size(images, 2) - 1)) / 10; 

%% If patient is Head First (HF) or Feet First (FF)

% Assume HFS
% info.ImageOrientationPatient = [1;0;0;0;1;0];

if isequal(info.ImageOrientationPatient, [1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [-1;0;0;0;-1;0]) 
    
    % Store patient position
    if info.ImageOrientationPatient(5) == 1
        img.position = 'HFS';
    elseif info.ImageOrientationPatient(5) == -1
        img.position = 'HFP';
    end

    % Sort sliceLocations vector in descending order
    [~, indices] = sort(sliceLocations, 'descend');
    
    % Store start voxel IEC-Y coordinate (in cm)
    img.start(3) = -max(sliceLocations) / 10;
    
elseif isequal(info.ImageOrientationPatient, [-1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [1;0;0;0;-1;0]) 
    
    % Store patient position
    if info.ImageOrientationPatient(5) == 1
        img.position = 'FFS';
    elseif info.ImageOrientationPatient(5) == -1
        img.position = 'FFP';
    end
    
    % Sort sliceLocations vector in ascending order
    [~,indices] = sort(sliceLocations, 'ascend');

    % Store start voxel IEC-Y coordinate (in cm)
    img.start(3) = min(sliceLocations) / 10;
    
end 

%% If the position wasn't identified, error
if ~isfield(img, 'position')
    
    error('The patient position could not be determined');
    
end

%% Compute slice location differences
widths = diff(sliceLocations(indices));

% Verify that slice locations do not differ significantly (1%)
if abs(max(widths) - min(widths))/mean(widths) > 0.01
    error(['Slice positions differ by more than 1%, suggesting ', ...
        'variable slice spacing. This is not supported.']);
end

%% Retrieve X/Z voxel widths from DICOM header (in cm)
img.width(1) = info.PixelSpacing(1) / 10;
img.width(2) = info.PixelSpacing(2) / 10;

%% Store mean slice position difference as IEC-Y width (in cm)
img.width(3) = abs(mean(widths)) / 10;

%% Initialize daily image data array as single type
img.data = single(zeros(size(images,3), size(images,2), size(images,1)));

%% Loop through each slice
for i = 1:length(sliceLocations)
    
    % Set the image data based on the index value
    img.data(:,:,i) = ...
        single(rot90(permute(images(indices(i),:,:), [2 3 1])));
    
end

%% Remove values less than zero 
% (some DICOM images place voxels outside the field of view 
% to negative values)
img.data = max(img.data, 0);

%% Flip images in IEC-X direction
img.data = flip(img.data, 1);

%% Create dimensions structure field based on the daily image size
img.dimensions = size(img.data);

%% Update waitbar
if exist('progress', 'var') && ishandle(progress)
    waitbar(1.0, progress, 'Image loading completed');
end

% If an image was unsuccessfully loaded
if ~isfield(img, 'dimensions')
    
    error('DICOM image data could not be parsed');
    
end

%% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

%% Clear temporary variables
clear i images info sliceLocations indices progress widths;

catch err
    
    % Delete progress handle if it exists
    if exist('progress', 'var') && ishandle(progress)
        delete(progress); 
    end
    
    % Catch errors and rethrow: 
    rethrow(err);
    
end

end