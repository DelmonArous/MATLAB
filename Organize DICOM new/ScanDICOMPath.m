function array = ScanDICOMPath(path, varargin)

% The function recursively searches a provided path for DICOM data
% and returns a cell array of DICOM images, RT structure sets, RT plan, and
% RT dose files along with basic header information found within the
% directory. DICOM files must contain the following minimum tags:
% SOPInstanceUID, MediaStorageSOPClassUID, Modality, PatientName, and
% PatientID.
%
% This function will display a progress bar while it loads unless MATLAB
% was executed with the -nodisplay, -nodesktop, or -noFigureWindows flags
% or if a 'Progress' input option is set to false.
%
% The following variables are required for proper execution:
%   path: string containing the path to the DICOM files, cell array of
%       files, or path to a single file
%
% Upon successful completion, the function will return an n x 12 cell
% array, where n is the number of files returned and the columns correspond
% to the following values:
%   Column 1: string containing the file name
%   Column 2: string containing the full path to the file
%   Column 3: string containing the file modality ('CT', 'MR', 'RTDOSE',
%       'RTPLAN', or 'RTSTRUCT')
%   Column 4: string containing the DICOM instance UID
%   Column 5: string containing the patient name, starting with the last
%       name separated by commas
%   Column 6: string containing the patient's ID
%   Column 7: string containing the study date
%   Column 8: string containing the frame of reference UID
%   Column 9: string containing the study or referenced study UID
%   Column 10: if CT/MR/PT, a string containing series description.
%       If RTPLAN or RTDOSE with a corresponding RTPLAN in the list,
%       a string containing the plan name
%   Column 11: if RTDOSE, the dose type ('PLAN' or 'BEAM')
%   Column 12: if RTPLAN, a cell array of beam names, or if a BEAM RTDOSE,
%       the corresponding beam name
%   Column 13: if RTPLAN, a cell array of machine names, or if RTDOSE, the
%       corresponding machine name
%   Column 14: if RTPLAN, a cell array of beam energies, or if RTDOSE, the
%       corresponding energy
%
% Below are examples of how this function is used:
%
%   % Scan test_directory for DICOM files
%   list = ScanDICOMPath('../test_directory');
%
%   % Re-scan, hiding the progress bar
%   list = ScanDICOMPath('../test_directory', 'Progress', false);
%
%   % List all CT file names found within the list
%   list(ismember(t(:,3), 'CT'))
%
%   % List all unique Frame of Reference UIDs in the list
%   unique(list(:,8))

%% Set default options
opt.Progress = true;

%% Parse provided options
for i = 2:2:nargin
    opt.(varargin{i-1}) = varargin{i};
end

%% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows') && opt.Progress
    progress = waitbar(0, 'Scanning path for DICOM files');
end

%% Set list based on format of provided files
if iscell(path)
    list = path;
elseif isfolder(path)
    list = dir(fullfile(path, '**'));
else
    % If not, throw an error
    error('The provided input is not a folder');
end

%% Initialize return array of DICOM files
array = cell(0,14);

%% Loop through each folder, subfolder
for i = 1:length(list)
    
    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && opt.Progress
        waitbar(i/length(list), progress);
    end
    
    %% If the folder content is . or .., skip to next folder in list
    if strcmp(list(i).name, '.') || strcmp(list(i).name, '..')
        continue
    elseif list(i).isdir == 1
        % Otherwise, if the folder content is a subfolder
        continue;
    else % Otherwise, see if the file is a DICOM file
        
        %% Separate file path, name
        [filepath, filename, ext] = fileparts(fullfile(list(i).folder, list(i).name));
        
        %% Attempt to parse the DICOM header
        try
            % Execute dicominfo
            info = dicominfo(fullfile(list(i).folder, list(i).name));
            
            % Verify storage class field exists
            if ~isfield(info, 'MediaStorageSOPClassUID')
                continue
            end
            
            % Store basic contents
            new = horzcat([filename ext], filepath, info.Modality, ...
                info.SOPInstanceUID, ...
                strjoin(struct2cell(info.PatientName), ', '), ...
                info.PatientID, ...
                datestr(datenum(strrep(info.StudyDate,'.',''), ...
                'yyyymmdd'), 'yyyy-mm-dd'), cell(1,7));
            
            %% If CT, MR or PET, respectively
            if strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.2') || ...
                    strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.4') || ...
                    strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.128')
                
                % Verify that enhanced contents exist
                if isfield(info, 'FrameOfReferenceUID')
                    new{8} = info.FrameOfReferenceUID;
                end
                if isfield(info, 'StudyInstanceUID')
                    new{9} = info.StudyInstanceUID;
                end
                if isfield(info, 'SeriesDescription')
                    new{10} = regexprep(strrep(info.SeriesDescription,'*','star'), ...
                        '[\\/?<>|:"]', '-');
                end
                
                % Append to table array
                array = [array; new];
                
                %% Otherwise, if structure
            elseif strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.481.3')
                
                % Verify that enhanced contents exist
                if isfield(info, 'ReferencedFrameOfReferenceSequence')
                    new{8} = info.ReferencedFrameOfReferenceSequence.Item_1...
                        .FrameOfReferenceUID;
                end
                if isfield(info, 'ReferencedStudySequence')
                    new{9} = info.ReferencedStudySequence...
                        .Item_1.ReferencedSOPInstanceUID;
                end
                
                % Append to table array
                array = [array; new];
                
                %% Otherwise, if dose
            elseif strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.481.2')
                
                % Verify that enhanced contents exist
                if isfield(info, 'FrameOfReferenceUID')
                    new{8} = info.FrameOfReferenceUID;
                end
                if isfield(info, 'ReferencedStudySequence')
                    new{9} = info.ReferencedStudySequence.Item_1...
                        .ReferencedSOPInstanceUID;
                end
                if isfield(info, 'ReferencedRTPlanSequence')
                    new{10} = info.ReferencedRTPlanSequence...
                        .Item_1.ReferencedSOPInstanceUID;
                    if isfield(info.ReferencedRTPlanSequence...
                            .Item_1, 'ReferencedFractionGroupSequence')
                        new{12} = {info.ReferencedRTPlanSequence...
                            .Item_1.ReferencedFractionGroupSequence.Item_1...
                            .ReferencedFractionGroupNumber
                            info.ReferencedRTPlanSequence...
                            .Item_1.ReferencedFractionGroupSequence.Item_1...
                            .ReferencedBeamSequence.Item_1...
                            .ReferencedBeamNumber};
                    end
                end
                if isfield(info, 'DoseSummationType')
                    new{11} = info.DoseSummationType;
                end
                
                % Append to table array
                array = [array; new];
                
                %% Otherwise, if RT plan
            elseif strcmp(info.MediaStorageSOPClassUID, ...
                    '1.2.840.10008.5.1.4.1.1.481.5')
                
                % Verify that enhanced contents exist
                if isfield(info, 'FrameOfReferenceUID')
                    new{8} = info.FrameOfReferenceUID;
                end
                if isfield(info, 'ReferencedStudySequence')
                    new{9} = info.ReferencedStudySequence...
                        .Item_1.ReferencedSOPInstanceUID;
                end
                if isfield(info, 'RTPlanLabel')
                    new{10} = info.RTPlanLabel;
                end
                if isfield(info, 'FractionGroupSequence') && ...
                        isfield(info, 'BeamSequence')
                    
                    % Store referenced beam numbers as 2D cell array of
                    % fraction group and beam numbers
                    new{12} = cell(0);
                    groups = fieldnames(info.FractionGroupSequence);
                    for j = 1:length(groups)
                        beams = fieldnames(info.FractionGroupSequence...
                            .(groups{j}).ReferencedBeamSequence);
                        for k = 1:length(beams)
                            new{12}{info.FractionGroupSequence.(groups{j})...
                                .FractionGroupNumber, k} = info...
                                .FractionGroupSequence...
                                .(groups{j}).ReferencedBeamSequence...
                                .(beams{k}).ReferencedBeamNumber;
                        end
                    end
                    
                    % Replace beam numbers with names and add machine/energy
                    new{13} = cell(size(new{12}));
                    new{14} = cell(size(new{12}));
                    beams = fieldnames(info.BeamSequence);
                    for j = 1:length(beams)
                        for k = 1:size(new{12},1)
                            for l = 1:size(new{12},2)
                                if ( isnumeric(new{12}{k,l}) && ...
                                        new{12}{k,l} == info.BeamSequence...
                                        .(beams{j}).BeamNumber )
                                    new{12}{k,l} = info.BeamSequence...
                                        .(beams{j}).BeamName;
                                    new{13}{k,l} = info.BeamSequence...
                                        .(beams{j}).TreatmentMachineName;
                                    new{14}{k,l} = num2str(info.BeamSequence...
                                        .(beams{j}).ControlPointSequence...
                                        .Item_1.NominalBeamEnergy);
                                    
                                    % Add non-standard identifier
                                    if ( isfield(info.BeamSequence...
                                            .(beams{j}), 'FluenceMode') && ...
                                            isfield(info.BeamSequence...
                                            .(beams{j}), 'FluenceModeID')&& ...
                                            strcmp(info.BeamSequence.(beams{j})...
                                            .FluenceMode, 'NON_STANDARD') )
                                        new{14}{k,l} = [new{14}{k,l}, ...
                                            info.BeamSequence.(beams{j})...
                                            .FluenceModeID];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Append to table array
                array = [array; new];
            end
            
        catch
            % Otherwise, the file is either corrupt or not a real DICOM
            % file, so throw an error
            warning(['File ', filename, ' is not a valid DICOM object.' ...
                newline 'Directory: ', filepath]);
            continue;
            
        end
    end
end

%% Loop through list, matching RT DOSE volumes to RT PlANs
for i = 1:size(array,1)
    
    % Match plan UID to UID of other DICOM files
    if ~isempty(array{i,10})
        [m,idx] = ismember(array{i,10}, array(:,4));
        
        % If RTDOSE is BEAM, match beam name, machine, and energy
        if m && strcmp(array{i,3}, 'RTDOSE') && ...
                strcmp(array{i,11}, 'BEAM') && ~isempty(array{i,12})
            array{i,14} = array{idx,14}{array{i,12}{1}, array{i,12}{2}};
            array{i,13} = array{idx,13}{array{i,12}{1}, array{i,12}{2}};
            array{i,12} = array{idx,12}{array{i,12}{1}, array{i,12}{2}};
        end
        
        % If RTDOSE is PLAN, duplicate PLAN cell arrays
        if m && strcmp(array{i,3}, 'RTDOSE') && strcmp(array{i,11}, 'PLAN')
            array{i,14} = array{idx,14};
            array{i,13} = array{idx,13};
            array{i,12} = array{idx,12};
        end
        
        % If RTDOSE
        if m && strcmp(array{i,3}, 'RTDOSE')
            array{i,10} = array{idx,10};
        end
    end
end

%% Close waitbar
if exist('progress', 'var') == 1 && ishandle(progress)
    close(progress);
end

%% Clear temporary variables
clear i j k l filepath filename ext list new info progress;