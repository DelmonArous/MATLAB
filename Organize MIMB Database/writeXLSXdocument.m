function [] = writeXLSXdocument(xlsxpath, data)

[filepath, filename, ext] = fileparts(xlsxpath);

try
    % If xlsx writing is successful, store the information
    [~, ~, data] = xlsread(fullfile(filepath, [filename ext]));
catch
    % Otherwise, the file is either corrupt or not a real .csv
    % file, so throw an error
    warning(['File ', filename, ' is not a valid .xlsx file.' ...
        newline 'Directory: ', filepath])
end

end