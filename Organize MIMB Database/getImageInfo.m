function [imginfo] = getImageInfo(path)

% The function stores relevant DICOM info of the images in the directory
% 'path' and returns the information stored in the object structure
% 'imginfo'

% The following variables are required for proper execution:
%   path: string containing the path to the DICOM files

fileList = getAllFiles(path);
% path

for i = 1:numel(fileList)
    
    [filepath, filename, ~] = fileparts(fileList{i});
    imginfo{i}.FilePath = filepath;
    imginfo{i}.FileName = filename;
    %     fileList{i}
    %     filepath
    %     filename
    
    try
        % If dicominfo is successful, store the header information
        info = dicominfo(fullfile(filepath, filename), 'UseDictionaryVR', true);
        [~, MSGID] = lastwarn();
        warning('off', MSGID)
    catch err
        if any(err.identifier) % strcmp(exception.identifier, 'MATLAB:structRefFromNonStruct')
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
            'TEST22222'
            % Otherwise, the file is either corrupt or not a real DICOM
            % file, so throw an error
            warning(['File ', filename, ...
                ' is not a valid DICOM object. Directory: ', filepath]);
            
            %         fileID = fopen(fullfile(path, 'Pat_warn.txt'), 'a');
            %         fprintf(fileID, 'File %s is not a valid DICOM object.\nDirectory: %s\n', ...
            %             filename, filepath);
            %         fclose(fileID);
            continue;
        end
        
    end
    
    %% Read Unique Identifier (UID) and patient characteristics
    imginfo{i}.SOPClassUID = info.SOPClassUID;
    imginfo{i}.StudyInstanceUID = info.StudyInstanceUID;
    imginfo{i}.SeriesInstanceUID = info.SeriesInstanceUID;
    
    if isfield(info, 'StudyDate')
        imginfo{i}.StudyDate = ['StudyDate ' datestr(datenum(...
            strrep(info.StudyDate,'.',''),'yyyymmdd'), 'yyyy-mm-dd')];
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
        
        if strcmp(info.Modality, 'MR')
            imginfo{i}.Modality = 'MR (Magnetic Resonance)';
        elseif strcmp(info.Modality, 'PT')
            imginfo{i}.Modality = 'PET (Positron Emission Tomography)';
        elseif strcmp(info.Modality, 'AR')
            imginfo{i}.Modality = 'AR (Autorefraction)';
        elseif strcmp(info.Modality, 'ASMT')
            imginfo{i}.Modality = 'ASMT (Content Assessment Results)';
        elseif strcmp(info.Modality, 'AU')
            imginfo{i}.Modality = 'AU (Audio)';
        elseif strcmp(info.Modality, 'BDUS')
            imginfo{i}.Modality = 'BDUS (Bone Densitometry (ultrasound))';
        elseif strcmp(info.Modality, 'BI')
            imginfo{i}.Modality = 'BI (Biomagnetic Imaging)';
        elseif strcmp(info.Modality, 'BMD')
            imginfo{i}.Modality = 'BMD (Bone Densitometry (X-Ray))';
        elseif strcmp(info.Modality, 'CR')
            imginfo{i}.Modality = 'CR (Computed Radiography)';
        elseif strcmp(info.Modality, 'CT')
            imginfo{i}.Modality = 'CT (Computed Tomography)';
        elseif strcmp(info.Modality, 'DG')
            imginfo{i}.Modality = 'DG (Diaphanography)';
        elseif strcmp(info.Modality, 'DOC')
            imginfo{i}.Modality = 'DOC (Document)';
        elseif strcmp(info.Modality, 'DX')
            imginfo{i}.Modality = 'DX (Digital Radiography)';
        elseif strcmp(info.Modality, 'ECG')
            imginfo{i}.Modality = 'ECG (Electrocardiography)';
        elseif strcmp(info.Modality, 'EPS')
            imginfo{i}.Modality = 'EPS (Cardiac Electrophysiology)';
        elseif strcmp(info.Modality, 'ES')
            imginfo{i}.Modality = 'ES (Endoscopy)';
        elseif strcmp(info.Modality, 'FID')
            imginfo{i}.Modality = 'FID (Fiducials)';
        elseif strcmp(info.Modality, 'GM')
            imginfo{i}.Modality = 'GM (General Microscopy)';
        elseif strcmp(info.Modality, 'CR')
            imginfo{i}.Modality = 'CR (Computed Radiography)';
        elseif strcmp(info.Modality, 'HC')
            imginfo{i}.Modality = 'HC (Hard Copy)';
        elseif strcmp(info.Modality, 'HD')
            imginfo{i}.Modality = 'HD (Hemodynamic Waveform)';
        elseif strcmp(info.Modality, 'IO')
            imginfo{i}.Modality = 'IO (Intra-Oral Radiography)';
        elseif strcmp(info.Modality, 'IOL') 
            imginfo{i}.Modality = 'IOL (Intraocular Lens Data)';
        elseif strcmp(info.Modality, 'IVOCT')
            imginfo{i}.Modality = 'IVOCT (Intravascular Optical Coherence Tomography)';
        elseif strcmp(info.Modality, 'IVUS')
            imginfo{i}.Modality = 'IVUS (Intravascular Ultrasound)';
        elseif strcmp(info.Modality, 'KER')
            imginfo{i}.Modality = 'KER (Keratometry)';
        elseif strcmp(info.Modality, 'KO')
            imginfo{i}.Modality = 'KO (Key Object Selection)';
        elseif strcmp(info.Modality, 'LEN')
            imginfo{i}.Modality = 'LEN (Lensometry)';
        elseif strcmp(info.Modality, 'LS')
            imginfo{i}.Modality = 'LS (Laser Surface Scan)';
        elseif strcmp(info.Modality, 'MG')
            imginfo{i}.Modality = 'MG (Mammography)';
        elseif strcmp(info.Modality, 'NM')
            imginfo{i}.Modality = 'NM (Nuclear Medicine)';
        elseif strcmp(info.Modality, 'OAM')
            imginfo{i}.Modality = 'OAM (Ophthalmic Axial Measurements)';
        elseif strcmp(info.Modality, 'OCT')
            imginfo{i}.Modality = 'OCT (Optical Coherence Tomography (non-Ophthalmic))';
        elseif strcmp(info.Modality, 'OP')
            imginfo{i}.Modality = 'OP (Ophthalmic Photography)';
        elseif strcmp(info.Modality, 'OPM')
            imginfo{i}.Modality = 'OPM (Ophthalmic Mapping)';
        elseif strcmp(info.Modality, 'OPT')
            imginfo{i}.Modality = 'OPT (Ophthalmic Tomography)';
        elseif strcmp(info.Modality, 'OT')
            imginfo{i}.Modality = 'OT (Other)';
        elseif strcmp(info.Modality, 'OPV')
            imginfo{i}.Modality = 'OPV (Ophthalmic Visual Field)';
        elseif strcmp(info.Modality, 'OSS')
            imginfo{i}.Modality = 'OSS (Optical Surface Scan)';
        elseif strcmp(info.Modality, 'PLAN')
            imginfo{i}.Modality = 'PLAN	(Plan)';
        elseif strcmp(info.Modality, 'PR')
            imginfo{i}.Modality = 'PR (Presentation State)';
        elseif strcmp(info.Modality, 'PX')
            imginfo{i}.Modality = 'PX (Panoramic X-Ray)';
        elseif strcmp(info.Modality, 'REG')
            imginfo{i}.Modality = 'REG (Registration)';
        elseif strcmp(info.Modality, 'CR')
            imginfo{i}.Modality = 'CR (Computed Radiography)';
        elseif strcmp(info.Modality, 'RESP')
            imginfo{i}.Modality = 'RESP	(Respiratory Waveform)';
        elseif strcmp(info.Modality, 'RF')
            imginfo{i}.Modality = 'RF (Radio Fluoroscopy)';
        elseif strcmp(info.Modality, 'RG')
            imginfo{i}.Modality = 'RG (Radiographic Imaging (conventional film/screen))';
        elseif strcmp(info.Modality, 'RTDOSE')
            imginfo{i}.Modality = 'RTDOSE (Radiotherapy Dose)';
        elseif strcmp(info.Modality, 'RTIMAGE')
            imginfo{i}.Modality = 'RTIMAGE (Radiotherapy Image)';
        elseif strcmp(info.Modality, 'RTPLAN')
            imginfo{i}.Modality = 'RTPLAN (Radiotherapy Plan)';
        elseif strcmp(info.Modality, 'RTRECORD')
            imginfo{i}.Modality = 'RTRECORD	(RT Treatment Record)';
        elseif strcmp(info.Modality, 'RTSTRUCT')
            imginfo{i}.Modality = 'RTSTRUCT	(Radiotherapy Structure Set)';
        elseif strcmp(info.Modality, 'RWV')
            imginfo{i}.Modality = 'RWV (Real World Value Map)';
        elseif strcmp(info.Modality, 'SEG')
            imginfo{i}.Modality = 'SEG (Segmentation)';
        elseif strcmp(info.Modality, 'SM')
            imginfo{i}.Modality = 'SM (Slide Microscopy)';
        elseif strcmp(info.Modality, 'SMR')
            imginfo{i}.Modality = 'SMR (Stereometric Relationship)';
        elseif strcmp(info.Modality, 'SR')
            imginfo{i}.Modality = 'SR (SR Document)';
        elseif strcmp(info.Modality, 'SRF')
            imginfo{i}.Modality = 'SRF (Subjective Refraction)';
        elseif strcmp(info.Modality, 'STAIN')
            imginfo{i}.Modality = 'STAIN (Automated Slide Stainer)';
        elseif strcmp(info.Modality, 'TG')
            imginfo{i}.Modality = 'TG (Thermography)';
        elseif strcmp(info.Modality, 'US')
            imginfo{i}.Modality = 'US (Ultrasound)';
        elseif strcmp(info.Modality, 'VA')
            imginfo{i}.Modality = 'VA (Visual Acuity)';
        elseif strcmp(info.Modality, 'XA')
            imginfo{i}.Modality = 'XA (X-Ray Angiography)';
        elseif strcmp(info.Modality, 'XC')
            imginfo{i}.Modality = 'XC (External-camera Photography)';
        else
            imginfo{i}.Modality = info.Modality;
        end
        
    end
    
end

clear fileList info filepath filename MSGID

end