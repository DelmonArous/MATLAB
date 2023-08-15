function [] = generateT2starMap(sourcepath, imgseq)

folderList = getAllFolders(sourcepath);

for i = 1:length(folderList)
    
    [folderpath, foldername, ~] = fileparts(folderList{i});
    
    if contains(foldername, 'FUNCPROST') && ...
            exist(fullfile(folderpath, [foldername '\' imgseq]), 'dir')
        spoiledGREfit_T2star(fullfile(folderpath, [foldername '\' imgseq]));
    end
    
end

clear folderList folderpath foldername i

end