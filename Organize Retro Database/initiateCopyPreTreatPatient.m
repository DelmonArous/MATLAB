function [] = initiateCopyPreTreatPatient(sourcepath, destpath, TreatDate)

[~, patID, ~] = fileparts(sourcepath);
fileList = getAllFiles(sourcepath);

flag = 'false';
strstruct = {''};
counter = 0;

for i = 1:numel(fileList)
    
    [filepath, filename, ~] = fileparts(fileList{i});
    
    try
        % If dicominfo is successful, store the header information
        info = dicominfo(fullfile(filepath, filename), 'UseDictionaryVR', true);
        [~, MSGID] = lastwarn();
        warning('off', MSGID)
    catch
        exception
        if any(exception.identifier) % strcmp(exception.identifier, 'MATLAB:structRefFromNonStruct')
            try
                info = dicominfo(fullfile(filepath, filename), 'UseDictionaryVR', false);
                [~, MSGID] = lastwarn();
                warning('off', MSGID)
            catch
                warning(['File ', filename, ...
                    ' is not a valid DICOM object. Directory: ', filepath]);
                continue;
            end
        else
            % Otherwise, the file is either corrupt or not a real DICOM
            % file, so throw an error
            warning(['File ', filename, ...
                ' is not a valid DICOM object. Directory: ', filepath]);
            
%             fileID = fopen(fullfile('C:\Users\Delmon\Desktop\Warning', ...
%                 strcat(patID, '.txt'), 'a'));
%             fprintf(fileID, 'File %s is not a valid DICOM object.\nDirectory: %s\n', ...
%                 filename, filepath);
%             fclose(fileID);
            continue;
        end
    end
    
    %% Store DICOM metadata for the study date
    % må returnere alle parameterene for denne pasienten
    % dette for å lage en array av alle pasientene
    % problem: flere sekvenser
    
    if isfield(info, 'StudyDate')
        StudyDate = datestr(datenum(strrep(info.StudyDate,'.',''),'yyyymmdd'), ...
            'yyyy-mm-dd');
    else
        warning(['File ', filename, ' has no study date.' ...
            newline 'Directory: ', filepath]);
        continue
    end
    
    %% Convert dates into datetime objects
    StudyDate = datetime(StudyDate, 'InputFormat', 'yyyy-MM-dd');
    TreatDate = datetime(TreatDate, 'InputFormat', 'yyyy-MM-dd');
    
    %% Get all pretreatment DICOM images
    if (StudyDate < TreatDate)
        
        flag = 'true';
        StudyDatestr = datestr(StudyDate);
        TreatDatestr = datestr(TreatDate);
        %newStudyDatestr = datestr(datenum(StudyDatestr, 'yyyy-mmm-dd'))
        %newTreatDatestr = datestr(datenum(TreatDatestr, 'yyyy-mmm-dd'))
        
        %% Store DICOM metadata from the header
        if isfield(info, 'SeriesDescription')
            SeriesDescription = regexprep(...
                strrep(info.SeriesDescription,'*','star'), '[\\/?<>|:"]', '-');
        end
        if isfield(info, 'ScanningSequence')
            ScanningSequence = info.ScanningSequence;
        end
        if isfield(info, 'MRAcquisitionType')
            MRAcquisitionType = info.MRAcquisitionType;
        end
        if isfield(info, 'SliceThickness')
            SliceThickness = info.SliceThickness;
        end
        if isfield(info, 'RepetitionTime')
            RepetitionTime = info.RepetitionTime;
        end
        if isfield(info, 'EchoTime')
            EchoTime = info.EchoTime;
        end
        if isfield(info, 'FlipAngle')
            FlipAngle = info.FlipAngle;
        end
        if isfield(info, 'Rows')
            Rows = info.Rows;
        end
        if isfield(info, 'Columns')
            Columns = info.Columns;
        end
        if isfield(info, 'PixelSpacing')
            PixelSpacing = info.PixelSpacing;
        end
        if isfield(info, 'Manufacturer')
            Manufacturer = info.Manufacturer;
        end
        if isfield(info, 'ManufacturerModelName')
            ManufacturerModelName = info.ManufacturerModelName;
        end
        
        %% Append header metadata to .txt file
        if ~any(strcmp(strstruct, SeriesDescription))
            
            strstruct{end+1} = SeriesDescription;
            counter = counter + 1;
            
            if exist(fullfile(destpath, 'info2.txt'), 'file') ~= 2
                fileID = fopen(fullfile(destpath, 'info2.txt'), 'w');
                fprintf(fileID, 'Patnr_MedInsightID,TreatDate,StudyDate,SeriesDescription,ScanningSequence,MRAcquisitionType,SliceThickness,TR,TE,FlipAngle,Matrix,PixelSpacing,Manufacturer,ManufacturerModelName\n');
            else
                fileID = fopen(fullfile(destpath, 'info2.txt'), 'a');
                
            end
            
            fprintf(fileID, '%s, %s, %s, %s, %s, %s, %f, %f, %f, %f, %dx%d, %f/%f, %s, %s\n', ...
                patID, TreatDatestr, StudyDatestr, SeriesDescription, ...
                ScanningSequence, MRAcquisitionType, SliceThickness, ...
                RepetitionTime, EchoTime, FlipAngle, Rows, Columns, ...
                PixelSpacing(1), PixelSpacing(2), ...
                Manufacturer, ManufacturerModelName);
            fclose(fileID);
            
        end
        %
        %         %% Create directory if the folder does not exist
        %         currentpath = fullfile(destpath, patID);
        %         if ~exist(currentpath, 'dir')
        %             mkdir(currentpath);
        %         end
        %
        %         currentpath = fullfile(currentpath, ...
        %             strcat('Pretreatment (TreatStart ', TreatDatestr, ')'));
        %         if ~exist(currentpath, 'dir')
        %             mkdir(currentpath);
        %         end
        %
        %         currentpath = fullfile(currentpath, ...
        %             strcat(SeriesDescription, ' (', StudyDatestr, ')')) ;
        %         if ~exist(fullfile(currentpath), 'dir')
        %             mkdir(currentpath);
        %         end
        %
        %         %% Copy file
        %         if ~strcmp(fullfile(filepath, filename), ...
        %                 fullfile(currentpath, filename))
        %             copyfile(fullfile(filepath, filename), ...
        %                 fullfile(currentpath, filename))
        %         end
        %
        %         % Reset destination
        %         currentpath = '';
        
    end
    
end

if strcmp(flag, 'false')
    fileID = fopen(fullfile('E:\Warning', strcat(patID, '.txt')), 'w');
    fprintf(fileID, '%s has not been imaged prior to treatment date %s', ...
        patID, datestr(TreatDate));
    fclose(fileID);
end

clear patID fileList flag filepath filename info fileID StudyDate ...
    TreatDate SeriesDescription currentpath

end