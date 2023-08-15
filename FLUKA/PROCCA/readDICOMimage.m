function [img] = readDICOMimage(path)

%% Get all files in the path directory
fileList = getAllFiles(path);
[~, name, ~] = fileparts(path);

%% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    % Start waitbar
    progress = waitbar(0, 'Loading DICOM images');
end

%% Initialize variables and read DICOM files

img.SOPClassUID             = '';
img.StudyInstanceUID        = '';
img.SeriesInstanceUID       = '';
img.FrameOfReferenceUID     = '';
img.MediaStorageSOPClassUID = '';
img.SOPInstanceUID          = cell(0);
img.PatientName             = '';
img.PatientID               = '';
img.PatientBirthDate        = '';
img.PatientSex              = '';
img.PatientAge              = '';
img.SeriesDescription       = '';
img.StudyDescription        = '';
info.InstitutionName        = '';

slicepositions              = [];
images                      = [];

for i = 1:numel(fileList)
    
    % Update progress bar
    if exist('progress', 'var') && ishandle(progress)
        waitbar(i/(numel(fileList) + 2), progress, ...
            sprintf('%s: Loading DICOM images (%i/%i)', ...
            name, i, numel(fileList)));
    end
    
    [filepath, filename, ext] = fileparts(fileList{i});
    
    if strcmp(ext, '.dcm')
        
        try
            % If dicominfo is successful, store the header information
            info = dicominfo(fullfile(filepath, filename), ...
                'UseDictionaryVR', true);
            [~, MSGID] = lastwarn();
            warning('off', MSGID)
        catch err
            
            if any(err.identifier)
                try
                    info = dicominfo(fullfile(filepath, filename), ...
                        'UseDictionaryVR', false);
                    [~, MSGID] = lastwarn();
                    warning('off', MSGID)
                catch
                    continue;
                end
            else
                % Otherwise, the file is either corrupt or not a real DICOM
                % file, so throw an error
                warning(['File ', filename, ' is not a valid DICOM object.' ...
                    newline 'Directory: ', filepath]);
                continue;
            end
            
        end
        
        % If this is the first DICOM image (and the class UID have not
        % yet been set
        if strcmp(img.SOPClassUID, '')
            
            % Store the UIDs, patient demographics
            img.SOPClassUID             = info.SOPClassUID;
            img.StudyInstanceUID        = info.StudyInstanceUID;
            img.SeriesInstanceUID       = info.SeriesInstanceUID;
            img.FrameOfReferenceUID     = info.FrameOfReferenceUID;
            
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
            if isfield(info, 'SeriesDescription')
                img.SeriesDescription = info.SeriesDescription;
            end
            if isfield(info, 'StudyDescription')
                img.StudyDescription = info.StudyDescription;
            end
            if isfield(info, 'InstitutionName')
                img.InstitutionName = info.InstitutionName;
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
        
        % Append this slice's image data to the images array
        if ~isempty(double(dicomread(info)))
            images(length(slicepositions) + 1,:,:) = dicomread(info);
            % double(dicomread(info));
        else
            warning(['Image array in file ', filename, ' is empty.' ...
                newline 'Directory: ', filepath]);
            continue;
        end
        
        % Add this slice's location to the slicepositions vector
        slicepositions(length(slicepositions) + 1) = ...
            info.ImagePositionPatient(3);
        
        % Add this slice's referenced image instance UIDs
        img.SOPInstanceUID{length(img.SOPInstanceUID) + 1} = ...
            info.SOPInstanceUID;
        
    end
    
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar((numel(fileList) + 1)/(numel(fileList) + 2), progress, ...
        sprintf('%s: Processing images', name));
end

%% Header modality tag
img.Modality                = info.Modality;
img.MediaStorageSOPClassUID = info.MediaStorageSOPClassUID;

%% Store start voxel IEC-X, IEC-Y and IEC-Z coordinates from DICOM header

% Retrieve start voxel IEC-X coordinate, in cm
img.start(1) = info.ImagePositionPatient(1) / 10 * ...
    info.ImageOrientationPatient(1);

% Adjust IEC-Z to inverted value, in cm
img.start(2) = -(info.ImagePositionPatient(2) * ...
    info.ImageOrientationPatient(5) + info.PixelSpacing(2) * ...
    (size(images, 2) - 1)) / 10;

% Retrieve X/Z voxel widths from DICOM header, in cm
img.width(1) = info.PixelSpacing(1) / 10;
img.width(2) = info.PixelSpacing(2) / 10;

%% Store image resolution and slicethickness from DICOM header
img.Rows            = info.Rows;
img.Columns         = info.Columns;
img.SliceThickness  = info.SliceThickness;

%% Set intercept and slope based on modality
if strcmp(info.MediaStorageSOPClassUID, '1.2.840.10008.5.1.4.1.1.2') && ...
        isfield(info, 'RescaleIntercept') && isfield(info, 'RescaleSlope')
%     img.RescaleIntercept   = info.RescaleIntercept;
%     img.RescaleSlope       = info.RescaleSlope;
    
    img.RescaleIntercept   = -1024;
    img.RescaleSlope       = 1;
elseif strcmp(info.MediaStorageSOPClassUID, '1.2.840.10008.5.1.4.1.1.2') && ...
        ~isfield(info, 'RescaleIntercept') && ~isfield(info, 'RescaleSlope')
    img.RescaleIntercept   = -1024;
    img.RescaleSlope       = 1;
else
    img.RescaleIntercept   = 0;
    img.RescaleSlope       = 1;
end

%% Window settings for CT
if strcmp(info.MediaStorageSOPClassUID, '1.2.840.10008.5.1.4.1.1.2') && ...
        isfield(info, 'WindowCenter') && isfield(info, 'WindowWidth')
   img.WindowCenter = info.WindowCenter;
   img.WindowWidth  = info.WindowWidth;
end
if isfield(info, 'WindowCenterWidthExplanation')
    img.WindowCenterWidthExplanation = info.WindowCenterWidthExplanation;
end

%% If patient orientation does not exist, assume HFS
if ~isfield(info, 'ImageOrientationPatient')
    info.ImageOrientationPatient = [1;0;0;0;1;0];
end

%% If patient is Head First
if isequal(info.ImageOrientationPatient, [1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [-1;0;0;0;-1;0])
    
    % Store patient position
    if info.ImageOrientationPatient(5) == 1
        img.PatientPosition = 'HFS';
        
    elseif info.ImageOrientationPatient(5) == -1
        img.PatientPosition = 'HFP';
    end
    
    % Sort sliceLocations vector in descending order
    [sortedslicepositions, indices] = sort(slicepositions, 'descend');
    
    % Store start voxel IEC-Y coordinate, in cm
    img.start(3) = -max(slicepositions) / 10;
    
    % Otherwise, if the patient is Feet First
elseif isequal(info.ImageOrientationPatient, [-1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [1;0;0;0;-1;0])
    
    % Store patient position
    if info.ImageOrientationPatient(5) == 1
        img.PatientPosition = 'FFS';
        
    elseif info.ImageOrientationPatient(5) == -1
        img.PatientPosition = 'FFP';
    end
    
    % % Sort sliceLocations vector in ascending order
    [sortedslicepositions, indices] = sort(slicepositions, 'ascend');
    
    % Store start voxel IEC-Y coordinate, in cm
    img.start(3) = min(slicepositions) / 10;
    
end

%% Compute slice location differences
widths = diff(slicepositions(indices));

% Verify that slice locations do not differ significantly (1%)
if abs(max(widths) - min(widths))/mean(widths) > 0.01
    error(['Slice positions differ by more than 1%, suggesting ', ...
        'variable slice spacing. This is not supported.']);
end

% Store mean slice position difference as IEC-Y width, in cm
img.width(3) = abs(mean(widths)) / 10;

%% Initialize daily image data array as single type
img.data = single(zeros(size(images, 3), size(images, 2), size(images, 1)));

%% Loop through each slice
for i = 1:length(slicepositions)
    % Set the image data based on the index value
    img.data(:, :, i) = ...
        double(rot90(permute(images(indices(i), :, :), [2 3 1])));
end

% Remove values less than zero (some DICOM images place voxels outside the
% field of view to negative values)
img.data = max(img.data, 0);

% Flip images in IEC-X direction
img.data = flip(img.data, 1);

% Create dimensions structure field based on the daily image size
img.dimensions = size(img.data);

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(1.0, progress, sprintf('%s: Image loading completed', name));
end

%% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

end