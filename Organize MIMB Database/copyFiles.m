function [] = copyFiles(patient, destpath)

% The function copies all files in the object 'patient' into a destination
% directory 'destpath'

% The following variables are required for proper execution:
%   patient: loadPatients object
%   destpath: string containing the destination path to where the DICOM
%               files will be copied

%% Create directory if the folder does not exist
if ~exist(fullfile(destpath), 'dir')
    mkdir(destpath);
end

%% Loop over patients and copy to destination directory
for i = 1:length(patient)
    
    currentpath = [destpath '\' patient{i}.name];
    if ~exist(fullfile(currentpath), 'dir')
        mkdir(currentpath);
    end
    
    if isfield(patient{i}, 'img')
        
        for j = 1:length(patient{i}.img)
            
            if isfield(patient{i}.img{j}, 'SeriesDescription')
                currentpath = [currentpath '\' patient{i}.img{j}.SeriesDescription ...
                    ' (' patient{i}.img{j}.StudyDate ')'];
                if ~exist(fullfile(currentpath), 'dir')
                    mkdir(currentpath);
                end
            end
            
            % Copy file to final destination directory
            if ~strcmp(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
                    fullfile(currentpath, patient{i}.img{j}.FileName))
                copyfile(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
                    fullfile(currentpath, patient{i}.img{j}.FileName))
            end
            
            % Reset destination
            currentpath = [destpath '\' patient{i}.name];
            
        end
        
    end
    
    % Delete empty subfolders
    flag = 1;
    while (flag)
        flag = deleteEmptyDir([destpath '\' patient{i}.name]);
    end
    
end

clear currentpath

end