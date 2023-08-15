function [ROIpatients] = readROIpatients(path)

fileList = getAllFiles(path);

counter = 0;
for i = 1:numel(fileList)
    
    [filepath, filename, ext] = fileparts(fileList{i});
    
    if strcmp(ext, '.roi')
        counter = counter + 1;
        
        [~, patientname, ~] = fileparts(filepath);
        ROIpatients{counter}.Patientname = patientname;
        ROIpatients{counter}.Filepath = filepath;
        ROIpatients{counter}.Filename = filename;
        
%         try
%             % If xlsx reading is successful, store the information
%             [~, ~, data] = xlsread(fullfile(filepath, [filename ext]));
%         catch
%             % Otherwise, the file is either corrupt or not a real .csv
%             % file, so throw an error
%             warning(['File ', filename, ' is not a valid .xlsx file.' ...
%                 newline 'Directory: ', filepath]);
%             continue;
%         end
%         
%         ROIpatients{counter}.ROIcenterX = cell2mat(data(2,9));
%         ROIpatients{counter}.ROIcenterY = cell2mat(data(2,10));
%         ROIpatients{counter}.SOPInstanceUID = cell2mat(data(2,17));
%         ROIpatients{counter}.SeriesInstanceUID = cell2mat(data(2,18));
%         ROIpatients{counter}.StudyInstanceUID = cell2mat(data(2,19));
%         ROIpatients{counter}.NumOfPoints = cell2mat(data(2,20));
%         ROIpatients{counter}.slicelocation = cell2mat(data(2,23));
%         
%         pxX_start = 24; pxY_start = 25;
%         mmX_start = 21; mmY_start = 22;
%         step = 5;
%         
%         ROIpatients{counter}.pxX = cell2mat(data(2, pxX_start:step: ...
%             (pxX_start + step*(ROIpatients{counter}.NumOfPoints - 1))));
%         ROIpatients{counter}.pxY = cell2mat(data(2, pxY_start:step: ...
%             (pxY_start + step*(ROIpatients{counter}.NumOfPoints - 1))));
%         ROIpatients{counter}.mmX = cell2mat(data(2, mmX_start:step: ...
%             (mmX_start + step*(ROIpatients{counter}.NumOfPoints - 1))));
%         ROIpatients{counter}.mmY = cell2mat(data(2, mmY_start:step: ...
%             (mmY_start + step*(ROIpatients{counter}.NumOfPoints - 1))));
      
                try
                    % If reading is successful, store the information
                    fileID = fopen(fullfile(filepath, [filename ext]), 'r');
                    data = textscan(fileID, '%f %f', 'Delimiter', ';');
                catch
                    % Otherwise, the file is either corrupt or not a real .csv
                    % file, so throw an error
                    warning(['File ', filename, ' is not a valid .roi file.' ...
                        newline 'Directory: ', filepath]);
                    continue;
                end
        
                ROIpatients{counter}.pxX = data{1}./3;
                ROIpatients{counter}.pxY = data{2}./3;
        
    end
    
end

clear fileList filepath filename fileID ext counter

end