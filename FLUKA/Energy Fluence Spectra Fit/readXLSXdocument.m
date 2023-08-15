function [data] = readXLSXdocument(xlsxpath)

[filepath, filename, ext] = fileparts(xlsxpath);

try
    % If .xlsx reading is successful, store the information
    [~, ~, data] = xlsread(fullfile(filepath, [filename ext]));
catch
    % Otherwise, the file is either corrupt or not a real .xlsx file,
    % so throw error
    warning(['File ', filename, ' is not a valid .xlsx file.' ...
        newline 'Directory: ', filepath])
end

end