function [flag] = deleteEmptyDir(path)

folderList = getAllFolders(path);
flag = 0;

for i = 1:numel(folderList)
    
    [folderpath, foldername, ~] = fileparts(folderList{i});
    
    if length(dir(folderList{i})) == 2
        rmdir(fullfile(folderpath, foldername));
        flag = 1;
    end
    
end

end