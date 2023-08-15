function [img] = readDICOMMR(path)

%% Get all files in the path directory
fileList = getAllFiles(path);
[~, name, ~] = fileparts(path);

%% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')  
    % Start progress bar
    progress = waitbar(0, sprintf('%s: Loading DICOM images', name)); 
end

%% Initialization
slicepositions = [];
imagedata = [];
TEvec = [];
metadata = {};

%% Read
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
            info = dicominfo(fullfile(filepath, filename), 'UseDictionaryVR', true);
        catch
            % Otherwise, the file is either corrupt or not a real DICOM
            % file, so throw an error
            warning(['File ', filename, ' is not a valid DICOM object.' ...
                newline 'Directory: ', filepath]);
            continue;
        end
        
        % Append this slice's image data to the images array
        if ~isempty(double(dicomread(info)))
            imagedata(:,:,length(slicepositions) + 1) = double(dicomread(info));
        else
            warning(['Image array in file ', filename, ' is empty.' ...
                newline 'Directory: ', filepath]);
            continue;
        end
        
        % Add this slice's metadata to the DICOM metadata vector
        metadata{length(metadata) + 1} = info;
        
        % Add this slice's echo time to the TEvec vector
        TEvec(length(TEvec) + 1) = info.EchoTime;
        
        % Add this slice's location to the slicepositions vector
        slicepositions(length(slicepositions) + 1) = info.ImagePositionPatient(3);
        
    end
    
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar((numel(fileList) + 1)/(numel(fileList) + 2), progress, ...
        sprintf('%s: Processing images', name));
end

%% Read series and study instance UID for the current patient
img.StudyInstanceUID = info.StudyInstanceUID;
img.SeriesInstanceUID = info.SeriesInstanceUID;

%% Read repetition time and flip angle of the MRI sequence
img.TR = info.RepetitionTime;
img.FA = info.FlipAngle;

%% Read start coordinates from DICOM header
img.ImagePositionPatient = info.ImagePositionPatient;

%% Read image resolution, pixel widths and slicethickness from DICOM header
img.Rows = info.Rows;
img.Columns = info.Columns;
img.PixelSpacing = info.PixelSpacing;
img.SliceThickness = info.SliceThickness;

%% If patient orientation does not exist, assume HFS
if ~isfield(info, 'ImageOrientationPatient')
    info.ImageOrientationPatient = [1;0;0;0;1;0];
end

info.ImageOrientationPatient = round(info.ImageOrientationPatient, 1);

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
    %[sortedslicepositions, indices] = sort(slicepositions, 'descend');
    [~, ~, indunique] = unique(slicepositions);
    %iterstep = 1;
    
elseif isequal(info.ImageOrientationPatient, [-1;0;0;0;1;0]) || ...
        isequal(info.ImageOrientationPatient, [1;0;0;0;-1;0])
    
    % Store patient position
    if (info.ImageOrientationPatient(5) == 1)
        img.position = 'FFS';
    elseif (info.ImageOrientationPatient(5) == -1)
        img.position = 'FFP';
    end
    
    % sliceposition vector sorted ascending
    %[sortedslicepositions, indices] = sort(slicepositions, 'ascend');
    [~, ~, indunique] = unique(slicepositions, 'stable');
    %iterstep = -1;
    
end

img.position

%% Number of slice locations (or unique slices)
img.N_slices = length(unique(slicepositions));

%% Loop through slices and sort based on index value
% slicepositions
% indunique
% a = find(indunique == 1)
% a = find(indunique == 2)
% a = find(indunique == 3)
% a = find(indunique == 4)
% 
% slicepositions
% TEvec
% size(imagedata)

for i = 1:img.N_slices % (img.N_slices - 1)
    
    img.slice{i}.TE = [];
    img.slice{i}.slicelocation = [];
    indices = find(indunique == i);
    
    if ~isempty(indices)
        
        img.slice{i}.info = metadata{indices(1)};
                
        for j = 1:length(indices) % indunique(i):iterstep:(indunique(i+1) - 1)
                   
            img.slice{i}.data(:,:,j) = imagedata(:,:,indices(j));
            img.slice{i}.TE = [img.slice{i}.TE; TEvec(indices(j))];
            img.slice{i}.slicelocation = [img.slice{i}.slicelocation; slicepositions(indices(j))];
                
        end
        
        [~, ind] = sort(img.slice{i}.TE);
        img.slice{i}.TE = img.slice{i}.TE(ind);
        img.slice{i}.data = img.slice{i}.data(:,:,ind);
        img.slice{i}.slicelocation = img.slice{i}.slicelocation(ind);
        
    end
        
end 

%% Remove image data values less than zero
% (some DICOM images place voxels outside the field of view to negative values)
for i = 1:img.N_slices
    img.slice{i}.data = max(img.slice{i}.data, 0);
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(1.0, progress, sprintf('%s: Image loading completed', name));
end

%% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

%% Clear temporary variables
clear i j progress N_slices iterstep TEvec imagedata slicepositions indices

end