function [] = organizeDICOM(patient, destpath)

% The function organizes the DICOM images allocated in the loadPatients 
% object 'patients' after patient-ID, modality, image creation date and image 
% series

% The following variables are required for proper execution:
%   patients: loadPatients object of relevant patients 

if ~exist(fullfile(destpath), 'dir')
    mkdir(destpath);
end

for i = 1:length(patient)
    
    currentpath = fullfile(destpath, patient{i}.name);
    currentpath
    if ~exist(fullfile(currentpath), 'dir')
        mkdir(currentpath);
    end
    
    for j = 1:length(patient{i}.img)
        
        if isfield(patient{i}.img{j}, 'Modality')
            currentpath = fullfile(currentpath, patient{i}.img{j}.Modality);
        else
            currentpath = fullfile(currentpath, 'Unknown modality');
        end
        if ~exist(fullfile(currentpath), 'dir')
            mkdir(currentpath);
        end
        
        if isfield(patient{i}.img{j}, 'StudyDate')
            currentpath = fullfile(currentpath, patient{i}.img{j}.StudyDate);
        else
            currentpath = fullfile(currentpath, 'Unknown study date');
        end
        if ~exist(fullfile(currentpath), 'dir')
            mkdir(currentpath);
        end
        
        if isfield(patient{i}.img{j}, 'SeriesDescription')
            currentpath = fullfile(currentpath, patient{i}.img{j}.SeriesDescription);
        else
            currentpath = fullfile(currentpath, 'Unknown series description');
        end
        if ~exist(fullfile(currentpath), 'dir')
            mkdir(currentpath);
        end
        
        % Move file to final destination directory
        % KOPIER HER ISTEDENFOR copyfile()!!
        if ~strcmp(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
                fullfile(currentpath, patient{i}.img{j}.FileName))            
            copyfile(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
                    fullfile(currentpath, patient{i}.img{j}.FileName))
%             movefile(fullfile(patient{i}.img{j}.FilePath, patient{i}.img{j}.FileName), ...
%                 fullfile(currentpath, patient{i}.img{j}.FileName));
        end
        
        % Reset destination
        currentpath = fullfile(destpath, patient{i}.name);
        
    end
    
    % Delete empty subfolders
%     flag = 1;
%     while (flag)
%         flag = deleteEmptyDir(patient{i}.dir);
%     end
    
end

clear currentpath flag

end