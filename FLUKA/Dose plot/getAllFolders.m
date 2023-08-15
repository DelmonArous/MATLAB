function [subfolders, filenames] = getAllFolders(path)

% The function recursively lists all subfolders and files under given folder 

% The following variables are required for proper execution: 
%   path: string containing the path to the DICOM files

% The following variables are returned upon succesful completion:
%   subfolders: stores all subfolders under given directory into a 
%               variable 'subfolders'
%   filenames: stores all filenames under given directory into a 
%               variable 'filenames'

if (nargin == 0)
    path = cd;
end

if (nargout == 1)
    subfolders = subfolder(path, '');
else
    [subfolders filenames] = subfolder(path, '', '');
end

    function [subfolders,filenames] = subfolder(path,subfolders,filenames)
        
        tmp = dir(path);
        tmp = tmp(~ismember({tmp.name},{'.' '..'}));
        
        for i = {tmp([tmp.isdir]).name}
            
            subfolders{end+1} = [path '\' i{:}];
            
            if (nargin == 2)
                subfolders = subfolder(subfolders{end},subfolders);
            else
                tmp = dir(subfolders{end});
                filenames{end+1} = {tmp(~[tmp.isdir]).name};
                [subfolders filenames] = subfolder(subfolders{end},subfolders,filenames);
            end
            
        end
    end
end