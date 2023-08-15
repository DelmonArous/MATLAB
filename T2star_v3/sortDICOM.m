function [] = sortDICOM(sourcepath, ROIpatients)

folderList = getAllFolders(sourcepath);

for i = 1:length(ROIpatients)
    
    ROIpatients{i}.Patientname
    currentpath = [sourcepath '\' ROIpatients{i}.Patientname '\ADC_new'];
    if ~exist(currentpath, 'dir')
        mkdir(currentpath);
    end
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        ROIpatients{i}.Patientname '\ADC']));
   
    if ~isempty(ind)
        
        img = readDICOMMR(folderList{ind});
        
        for j = 1:img.N_slices
            
            dicomwrite(uint16(img.slice{j}.data(:,:,1)), fullfile(currentpath, ...
                strcat('Image#', sprintf('%02d', j), '.dcm')), img.slice{j}.info, ...
                'CreateMode', 'Copy', 'CompressionMode', 'None', 'MultiframeSingleFile', 'true');
            
        end

    end
    
end

end