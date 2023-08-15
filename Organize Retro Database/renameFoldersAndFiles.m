function [] = renameFoldersAndFiles(path, n)

% The function allocates all folders (and subfolders) and files (and subfiles) 
% in the directory 'path'. Folders and files with names longer than
% 12 characters is renamed to the last n characters of the current folder/filername

% The following variables are required for proper execution:
%   path: string containing the path to the DICOM files
%   n: number of filename characters

fileList = getAllFiles(path);
folderList = getAllFolders(path);

counter = 0;
for i = 1:numel(fileList)
    
    % Get the directory, name and extension of the folder
    [filepath, filename, ext] = fileparts(fileList{i});
    
    if (length(filename) > 12) % (~isnan(str2double(ext)) || strcmpi(ext, '.dcm'))
 
        % New filename is the n last characters of filename
        newFilename = filename(end-(n-1):end);
        
        % If the newFilename exists, permute string, and then rename
        if exist(fullfile(filepath, newFilename), 'file') == 2
            
            % Find all possible permutations of the string newFilename
            permute = perms(newFilename);
            
            while (exist(fullfile(filepath, newFilename), 'file'))
                counter = counter + 1;
                if (counter <= length(permute))
                    newFilename = permute(counter, :);
                else
                    newFilename = num2str(round( ...
                        rand*(99999999 - 10000000) + 10000000));
                end
            end
            
            % Reset permutation index
            counter = 0;
            
        end
        
        % Rename the file and move to the directory filepath
        movefile(fullfile(filepath, [filename ext]), ...
            fullfile(filepath, newFilename));
        
    end
    
end

counter = 0;
for i = 1:numel(folderList)
    
    % Get the directory, name and extension of the folder
    [folderpath, foldername, ext] = fileparts(folderList{i});
    
    if (length(foldername) > 12)
        
        % New filename is the n last characters of filename
        newFoldername = foldername(end-(n-1):end);
        
        % If the newFoldername exists, permute string, and then rename
        if exist(fullfile(folderpath, newFoldername), 'file') == 7
            
            % Find all possible permutations of the string newFoldername
            permute = perms(newFoldername);
            
            while (exist(fullfile(folderpath, newFoldername), 'file'))
                counter = counter + 1;
                if (counter <= length(permute))
                    newFoldername = permute(counter, :);
                else
                    newFoldername = num2str(round( ...
                        rand*(99999999 - 10000000) + 10000000));
                end
            end
            
            % Reset permutation index
            counter = 0;
            
        end
        
        % Rename the folder and move to the same absolute path folderpath
        movefile(fullfile(folderpath, [foldername ext]), ...
            fullfile(folderpath, newFoldername));
        
    end
    
end

clear fileList folderList counter folderpath foldername filepath filename ...
    ext newFilename newFoldername

end