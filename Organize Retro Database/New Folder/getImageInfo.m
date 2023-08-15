function [imginfo] = getImageInfo(path)

% The function stores relevant DICOM info of the images in the directory 
% 'path' and returns the information stored in the object structure 
% 'imginfo'

% The following variables are required for proper execution:
%   path: string containing the path to the DICOM files

fileList = getAllFiles(path);

for i = 1:numel(fileList)
    
    [filepath, filename, ~] = fileparts(fileList{i});
    imginfo{i}.FilePath = filepath;
    imginfo{i}.FileName = filename;
    
    try
        % If dicominfo is successful, store the header information
        info = dicominfo(fullfile(filepath, filename), 'UseDictionaryVR', true);
    catch
        % Otherwise, the file is either corrupt or not a real DICOM
        % file, so throw an error
        warning(['File ', filename, ' is not a valid DICOM object.' ...
            newline 'Directory: ', filepath]);
        
%         fileID = fopen(fullfile(path, 'Pat_warn.txt'), 'a');
%         fprintf(fileID, 'File %s is not a valid DICOM object.\nDirectory: %s\n', ...
%             filename, filepath);
%         fclose(fileID);
        continue;
    end
    
    %% Read Unique Identifier (UID) and patient characteristics
    imginfo{i}.SOPClassUID = info.SOPClassUID;
    imginfo{i}.StudyInstanceUID = info.StudyInstanceUID;
    imginfo{i}.SeriesInstanceUID = info.SeriesInstanceUID;
   
    if isfield(info, 'StudyDate')
        imginfo{i}.StudyDate = datestr(datenum(...
            strrep(info.StudyDate,'.',''),'yyyymmdd'), 'yyyy-mm-dd');
    end
    if isfield(info, 'PatientName')
        imginfo{i}.PatientName = info.PatientName;
    end
    if isfield(info, 'PatientID')
        imginfo{i}.PatientID = info.PatientID;
    end
    if isfield(info, 'PatientBirthDate')
        imginfo{i}.PatientBirthDate = info.PatientBirthDate;
    end
    if isfield(info, 'PatientSex')
        imginfo{i}.PatientSex = info.PatientSex;
    end
    if isfield(info, 'PatientAge')
        imginfo{i}.PatientAge = info.PatientAge;
    end
    if isfield(info, 'StudyDescription')
        imginfo{i}.StudyDescription = info.StudyDescription;
    end
    if isfield(info, 'SeriesDescription')
        imginfo{i}.SeriesDescription = regexprep(...
            strrep(info.SeriesDescription,'*','star'), '[\\/?<>|:"]', '-');
    end
    
    %% Give object imginfo appropriate modality info
    if isfield(info, 'Modality')
        
        % MR = Magnetic Resonance
        % PT = Positron emission tomography (PET)
        % CT = Computed Tomography
        % CR = Computed Radiography
        % RTSTRUCT = Radiotherapy Structure Set
        % PR = Presentation State
        
        imginfo{i}.Modality = info.Modality;
      
        if strcmpi(info.Modality, 'MR')
            if isfield(info, 'ScanningSequence')
                imginfo{i}.ScanningSequence = info.ScanningSequence;
            end
            
            if isfield(info, 'MRAcquisitionType')
                imginfo{i}.MRAcquisitionType = info.MRAcquisitionType;
            end
        elseif strcmpi(info.Modality, 'PT')
            imginfo{i}.Modality = 'PET';
        end
        
    end
    
end

clear fileList info filepath filename 

end