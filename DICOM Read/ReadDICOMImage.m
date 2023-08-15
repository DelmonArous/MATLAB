function [img] = ReadDICOMImage(path, filenames)

% The function loads a series of single-frame DICOM CT or MR images and 
% returns a formatted structure for dose calculation.

% The following variables are required for proper execution: 
%   path: string containing the path to the DICOM files
%   filenames: cell array of strings containing all files to be loaded

% The following variables are returned upon succesful completion:
%   image: structure containing the image data, dimensions, width, type,
%       start coordinates, and key DICOM header values. The data is a three 
%       dimensional array of image values, while the dimensions, width, and 
%       start fields are three element vectors.  The DICOM header values 
%       are returned as a strings.

%% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    
    % Start progress bar
    progress = waitbar(0, 'Loading DICOM images');
    
end

%% Initialize empty variables for Unique Identifier (UID) and patient characteristics
img.SOPClassUID = '';
img.FrameOfReferenceUID = '';
img.StudyInstanceUID = '';
img.SeriesInstanceUID = '';
img.SOPInstanceUID = cell(0);
img.PatientName = '';
img.PatientID = '';
img.PatientBirthDate = '';
img.PatientAge = '';

%% Initialize empty 3D array for images and vector of slice locations
% (the data may not be loaded in correct order; these will be used to
% re-sort the slices later)
images = [];
slicepositions = [];

%% Loop through each file in filenames list
for i = 1:length(filenames)
    
    % Update progress bar
    if exist('progress', 'var') && ishandle(progress)
        waitbar(i/(length(filenames)+2), progress, ...
            sprintf('Loading DICOM images (%i/%i)', i, length(filenames)));
    end
    
    try 

        % Read the DICOM header of imagefile
        info = dicominfo(fullfile(path, filenames{i}));
    
    catch
        
        % Otherwise, warn that the file is either corrupt or not a real 
        % DICOM file, and then skip to next file in directory
        warning(['File ', filenames{i}, ' is not a valid DICOM image ', ...
            'and was skipped']);
        continue
        
    end
    
    % If this is the first DICOM image (and the class UID
    % have not yet been set), store UIDs and patient demographics
    if strcmp(img.SOPClassUID, '')
        
        img.SOPClassUID         = info.SOPClassUID;
        img.FrameOfReferenceUID = info.FrameOfReferenceUID;
        img.StudyInstanceUID    = info.StudyInstanceUID;
        img.SeriesInstanceUID   = info.SeriesInstanceUID;
        
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
        
    elseif ~strcmp(img.StudyInstanceUID, info.StudyInstanceUID)
        
        % Otherwise, if this file's study UID does not match the others,
        % multiple DICOM studies may be present in the same folder (not
        % currently supported)
        error(['Multiple DICOM Study Instance UIDs were found in ', ...
            'this list.  Please select only one study.']);
        
    elseif ~strcmp(img.SeriesInstanceUID, info.SeriesInstanceUID) 
        
        % Otherwise, if this file's series UID does not match the others,
        % multiple DICOM series may be present in the same folder (not
        % currently supported)
        error(['Multiple DICOM Series Instance UIDs were found in ', ...
            'this list.  Please select only one series.']);
        
    end
    
    % Add this slice's instance UID
    img.SOPInstanceUID{length(img.SOPInstanceUID) + 1} = info.SOPInstanceUID;
    
    % Add this slice's location to the slicepositions vector
    slicepositions(length(slicepositions)+1) = info.ImagePositionPatient(3);

    % Append this slice's image data to the images array
    images(size(images,1) + 1,:,:) = dicomread(info);
    
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar((length(filenames)+1)/(length(filenames)+2), progress, ...
        'Processing images');
end

%% Set image type based on DICOM header modality tag
img.type = info.Modality; 

%% Start IEC-X voxel coordinate from DICOM header (in cm)
img.start(1) = info.ImagePositionPatient(1) / 10 * ... 
    info.ImageOrientationPatient(1);

%% Start IEC-Z voxel coordinate (inverted value) (in cm)
img.start(2) = -(info.ImagePositionPatient(2) * ...
    info.ImageOrientationPatient(5) + info.PixelSpacing(2) * ...
    (size(images, 2) - 1)) / 10; 

%% If patient orientation does not exist, assume HFS
%if ~isfield(info, 'ImageOrientationPatient') 
    info.ImageOrientationPatient = [1;0;0;0;1;0];
%end

%% If the position wasn't identified, error
if ~isfield(img, 'position')
   
    error('The patient position could not be determined');
    
end

%% If patient is head first (HF) or feet first (FF)
if isequal(info.ImageOrientationPatient, [1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [-1;0;0;0;-1;0]) 
    
    % Store patient position
    if (info.ImageOrientationPatient(5) == 1)
        img.position = 'HFS';
    elseif (info.ImageOrientationPatient(5) == -1)
        img.position = 'HFP';
    end

    % slicepositions vector sorted in descending order
    [~, indices] = sort(slicepositions, 'descend');

    % Start voxel IEC-Y coordinate in cm
    img.start(3) = -max(slicepositions) / 10;

elseif isequal(info.ImageOrientationPatient, [-1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [1;0;0;0;-1;0]) 
    
    % Store patient position
    if (info.ImageOriefntationPatient(5) == 1)
        img.position = 'FFS';
    elseif (info.ImageOrientationPatient(5) == -1)
        img.position = 'FFP';
    end

    % sliceposition vector sorted ascending
    [~,indices] = sort(slicepositions, 'ascend');

    % Start voxel IEC-Y coordinate in cm
    img.start(3) = min(slicepositions) / 10;

end

%% Check if slice locations differ significantly (1%) 

% Estimate slice location difference
widthsdiff = diff(slicepositions(indices)); 
if abs(max(widthsdiff) - min(widthsdiff))/mean(widthsdiff) > 0.01
    error(['Error occured. \n', ...
        'Slice position vary more than 1% because of degree of variable slice spacing.']);
end

% Retrieve X/Z voxel widths from DICOM header, in cm
img.width(1) = info.PixelSpacing(1) / 10;
img.width(2) = info.PixelSpacing(2) / 10;
% Store mean slice position difference as IEC-Y width, in cm
img.width(3) = abs(mean(widthsdiff)) / 10;

%% Initialize daily image data array as single type
img.data = single(zeros(size(images, 3), size(images, 2), size(images, 1)));

%% Loop through slices and sort based on index value
for i = 1:length(slicepositions)
    img.data(:,:,i) = single(rot90(permute(images(indices(i),:,:), ...
        [2 3 1])));
end

%% Remove image data values less than zero 
% (some DICOM images place voxels outside the 
% field of view to negative values)
img.data = max(img.data, 0);

%% Flip images in IEC-X direction
img.data = flip(img.data, 1);

%% Create dimensions structure field based on the daily image size
img.dimensions = size(img.data);

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(1.0, progress, 'Image loading completed');
end

%% If an image was not successfully loaded
if ~isfield(image, 'dimensions')
    error('DICOM image data could not be parsed');
end

%% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

%% Clear temporary variables
clear i images info slicepositions indices progress widthsdiff;

end