function [fileList, folderList] = getAllFiles(path)

dirData = dir(path);                    % Get the data for the current directory
dirIndex = [dirData.isdir];             % Find the index for directories
fileList = {dirData(~dirIndex).name}';  % Get a list of the files

%% Prepend path to files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(path,x), fileList, 'UniformOutput', false);
end

%% Get a list of the subdirectories
subDirs = {dirData(dirIndex).name};

%% Find index of subdirectories that are not '.' or '..'
validIndex = ~ismember(subDirs, {'.', '..'});

%% Loop over valid subdirectories
for iDir = find(validIndex)
    
    % Get the subdirectory path
    nextDir = fullfile(path,subDirs{iDir});
    % Recursively call getAllFiles
    fileList = [fileList; getAllFiles(nextDir)];
    
end

end