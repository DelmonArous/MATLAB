function [struct] = ReadDICOMRTSS(varargin)

% This function reads a DICOM RT Structure Set (RTSS) file and
% extracts the structure information into a MATLAB cell array.
% This function will display a progress bar while it loads (unless MATLAB 
% was executed with the -nodisplay, -nodesktop, or -noFigureWindows flags).

% The following variables are required for proper execution: 
%   varargin{1}: string containing the path to the DICOM files
%   varargin{2}: string containing the DICOM RTSS file to be loaded
%   varargin{3}: structure of reference image.  Must include a frameRefUID
%       field referencing the structure set, as well as dimensions, width, 
%       and start fields. See LoadDICOMImages for more information.
%   varargin{4} (optional): flag indicating whether to ignore frame of
%       reference (1) or to verify it matches the CT (0, default)

% The following variable is returned upon succesful completion:
%   struct: cell array of structure names, color, start, width,
%       dimensions, frameRefUID, and 3D mask array of same size as
%       reference image containing fraction of voxel inclusion in structure


%% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    % Start progress bar
    progress = waitbar(0, 'Loading DICOM structures'); 
end

%% Read DICOM header info from file
info = dicominfo(fullfile(varargin{1}, varargin{2}));

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(0.1, progress);
end

%% Initialize structures
struct = cell(length(fieldnames(info.ROIContourSequence)) + 1, 1);

% Error if the DICOM file contains no structures
if isempty(struct)
    error('Error occured. No contours were found in the file.');
end

%% Loop through items in StructureSetROISequence
for item = fieldnames(info.StructureSetROISequence)'
    
    % Store number of ROI (contour) 
    ROINumber = info.StructureSetROISequence.(item{1}).ROINumber + 1;
    
    % Store name of ROI (contour)
    ROIName = info.StructureSetROISequence.(item{1}).ROIName;
    
    % If the structure frame of refererence matches the image 
    % frame of reference
    if strcmp(varargin{3}.FrameOfReferenceUID, info.StructureSetROISequence.(...
            item{1}).ReferencedFrameOfReferenceUID)
        
        % Store contour (ROI) name 
        struct{ROINumber}.ROIName = ROIName;
        
        % Store frame of reference UID for the contour
        struct{ROINumber}.frameRefUID = info.StructureSetROISequence.(...
            item{1}).ReferencedFrameOfReferenceUID;
        
    else
        
        warning(['The frame of reference for structure ', ROIName,...
            ' is not equal to that of the image.']);
        
    end
    
    % Clear temporary variables
    clear ROIName;
    
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(0.2, progress);
end

%% Administer image orientation (rotation vector) based on patient position
if isfield(varargin{3}, 'position')
    
    if strcmpi(varargin{3}.position, 'HFS')
        rotvec = [1,1,1];
    elseif strcmpi(varargin{3}.position, 'HFP')
        rotvec = [-1,-1,1];
    elseif strcmpi(varargin{3}.position, 'FFS')
        rotvec = [-1,1,-1];
    elseif strcmpi(varargin{3}.position, 'FFP')
        rotvec = [1,-1,-1];
    end

else
    
    % Assume HFS if patient position is indefined
    rotvec = [1,1,1];

end

%% Loop through items in ROIContourSequence

% Counter
counter = 0;

for item = fieldnames(info.ROIContourSequence)'
    
    % Increment
    counter = counter + 1;
    
    % Store number of contour (ROI)
    if isfield(info.ROIContourSequence.(item{1}), 'ReferencedROINumber')
        ReferencedROINumber = info.ROIContourSequence.(item{1}).ReferencedROINumber;
    else
        ReferencedROINumber = counter;
    end
    
    % Update progress bar
    if exist('progress', 'var') && ishandle(progress)
        waitbar(0.2 + 0.8*ReferencedROINumber/length(fieldnames(info.ROIContourSequence)), ...
            progress, sprintf('Loading DICOM structures (%i/%i)', ...
            ReferencedROINumber, length(fieldnames(info.ROIContourSequence))));
    end
    
    % If structure has a ROIName
    if isfield(struct{ReferencedROINumber}, 'ROIName')
        
        % Store the ROI color if available
        if isfield(info.ROIContourSequence.(item{1}), 'ROIDisplayColor')
            struct{ReferencedROINumber}.color = ...
                info.ROIContourSequence.(item{1}).ROIDisplayColor';
        else
            struct{ReferencedROINumber}.color = [0 0 0];
        end
        
        % Generate empty logical mask of the same size as the reference image
        struct{ReferencedROINumber}.mask = false(varargin{3}.dimensions);
        
        % Initialize structure volume and structure datapoint cell array
        struct{ReferencedROINumber}.volume = 0;
        struct{ReferencedROINumber}.datapoints = cell(0);
        
        % If a contour sequence does not exist, skip this structure
        if ~isfield(info.ROIContourSequence.(item{1}), 'ContourSequence')
            continue;
        end
        
        % Loop through each ContourSequence
        for cs = fieldnames(info.ROIContourSequence.(...
                item{1}).ContourSequence)'
            
            % Skip to next ContourSequence if no contour points exist,
            % otherwise load contour points
            if (info.ROIContourSequence.(item{1}).ContourSequence.(...
                    cs{1}).NumberOfContourPoints == 0)
                continue;
            else
                
                % Read contour datapoints and convert from mm to cm,
                % and reorder based on rotation vector
                datapoints = reshape(info.ROIContourSequence.(...
                    item{1}).ContourSequence.(cs{1}).ContourData, 3, [])' / 10;
                datapoints = datapoints.*repmat(rotvec, size(datapoints, 1), 1);
                
                % If points are empty, warn and continue
                if isempty(datapoints)
                    
                    warning(['Structure ', struct{ReferencedROINumber}.ROIName, ...
                        ' contains an empty contour']);
                    continue;
                    
                end
                
                % Store contour data
                struct{ReferencedROINumber}.datapoints{length(...
                    struct{ReferencedROINumber}.datapoints)+1} = datapoints;
                
                % Determine slice index by searching IEC-Y index using
                % nearest neighbor interpolation
                slice = interp1(varargin{3}.start(3):varargin{3}.width(3): ...
                    varargin{3}.start(3) + (varargin{3}.dimensions(3) - 1) ...
                    *varargin{3}.width(3), 1:varargin{3}.dimensions(3), ...
                    -datapoints(1,3), 'nearest', 0);
                
                % If the slice index is within the reference image
                if (slice ~= 0)
                    
                    % Convert a ROI polygon from datapoints to a binary 
                    % ROI mask. The voxels in mask that are outside and 
                    % inside the polygon are 0 and 1, respectively.
                    % Test if voxel centers are within polygon defined by 
                    % point data, adding result to structure mask.
                    mask = poly2mask((datapoints(:,2) + (varargin{3}.width(2)* ...
                        (varargin{3}.dimensions(2) - 1) + ...
                        varargin{3}.start(2)) + varargin{3}.width(2)/2) / ...
                        varargin{3}.width(2) + 1, ...
                        (datapoints(:,1) + varargin{3}.width(1)/2 - ...
                        varargin{3}.start(1)) / varargin{3}.width(1) + 1, ...
                        varargin{3}.dimensions(1), varargin{3}.dimensions(2));

                    % Voxels that are surrounded by even numbers of curves 
                    % are outside of the structure. 
                    % Subtract if the new mask overlaps with an existing 
                    % value, otherwise, add
                    if (max(max(mask + struct{ReferencedROINumber}.mask(:,:,slice))) == 2)
                        struct{ReferencedROINumber}.mask(:,:,slice) = ...
                            struct{ReferencedROINumber}.mask(:,:,slice) - mask;
                    else
                        struct{ReferencedROINumber}.mask(:,:,slice) = ...
                            struct{ReferencedROINumber}.mask(:,:,slice) + mask;
                    end
                    
                % Warn that the contour did not match a slice  
                else
                    
                    warning(['Structure ', struct{ReferencedROINumber}.ROIName, ...
                        ' contains contours outside of image array']);
                    
                end
                
            end
        end
        
        % Compute volume of the structure from the binary mask, however 
        % partial voxels are not included in the estimate
        struct{ReferencedROINumber}.volume = ...
            nnz(struct{ReferencedROINumber}.mask)*prod(varargin{3}.width);
        
        % Copy structure width, start, and dimensions arrays from image
        struct{ReferencedROINumber}.width = varargin{3}.width;
        struct{ReferencedROINumber}.start = varargin{3}.start;
        struct{ReferencedROINumber}.dimensions = varargin{3}.dimensions;
        
        % Warn if the mask is empty and clear structure
        if max(max(max(struct{ReferencedROINumber}.mask))) == 0
            
            warning(['Structure ', struct{ReferencedROINumber}.ROIName, ...
                ' is less than one voxel and will be removed']);
            struct{ReferencedROINumber} = [];
            
        end
        
    end
end

%% Remove empty cells in the structure
struct = struct(~cellfun('isempty', struct));

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(1.0, progress, 'Structure set loading completed');
end

%% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

%% Clear temporary files
clear info ReferencedROINumber counter item cs datapoints slice mask progress rotvec


%% Delete progress handle if it exists
if exist('progress', 'var') && ishandle(progress), delete(progress); end

   
end