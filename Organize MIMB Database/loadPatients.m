function [patient] = loadPatients(path)

% The function creates a structure of all the patients with DICOM images
% in the directory 'path'

% The following variables are required for proper execution:
%   path: string containing the path to the relevant patients

folderList = getAllFolders(path);

counter = 0;
for i = 1:length(folderList)
    
    [folderpath, foldername, ~] = fileparts(folderList{i});
    
    if strcmp(foldername, 'DICOM')
        counter = counter + 1;
        [path, ~, ~] = fileparts(folderpath);
        [~, name, ~] = fileparts(path);
        fullfile(folderpath, foldername)
        patient{counter}.name = name;
        patient{counter}.img = ...
            getImageInfo(fullfile(folderpath, foldername));
        patient{counter}.dir = name; % folderpath; % fullfile(folderpath, foldername);
    end
    
end

clear folderList folderpath foldername name ext i counter

end