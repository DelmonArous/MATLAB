function [] = organizeDICOM(patient)

% The function organizes the DICOM images allocated in the loadPatients 
% object 'patients' after patient-ID, modality, image creation date and image 
% series

% The following variables are required for proper execution:
%   patients: loadPatients object of relevant patients 

for i = 1:length(patient)
    
    currentpath = patient{i}.dir;
    if ~exist(fullfile(currentpath), 'dir')
        mkdir(currentpath);
    end
    
    for j = 1:length(patient{i}.img)
        
        currentpath = [currentpath '\' patient{i}.img{j}.Modality];
        if ~exist(fullfile(currentpath), 'dir')
            mkdir(currentpath);
        end
        
        currentpath = [currentpath '\' patient{i}.img{j}.StudyDate];
        if ~exist(fullfile(currentpath), 'dir')
            mkdir(currentpath);
        end
        
        if isfield(patient{i}.img{j}, 'SeriesDescription')
            currentpath = [currentpath '\' patient{i}.img{j}.SeriesDescription];
            if ~exist(fullfile(currentpath), 'dir')
                mkdir(currentpath);
            end
        end
        
        % Move file to final destination directory
        if ~strcmp(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
                fullfile(currentpath, patient{i}.img{j}.FileName))
            movefile(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
                fullfile(currentpath, patient{i}.img{j}.FileName));
        end
        
        % Reset destination
        currentpath = patient{i}.dir;
        
    end
    
    % Delete empty subfolders
    flag = 1;
    while (flag)
        flag = deleteEmptyDir(patient{i}.dir);
    end
    
end

clear currentpath flag

end