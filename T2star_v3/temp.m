function [] = temp(sourcepath, imgseq, ROIpatients)

folderList = getAllFolders(sourcepath);

%% Loop over each ROI for each patient
for i = 1:length(ROIpatients)
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        ROIpatients{i}.Patientname '\' imgseq]));
   
    if ~isempty(ind)
       
        fileList = getAllFiles(folderList{ind});
        
        for j = 1:numel(fileList)
            
            [p, n, e] = fileparts(fileList{j});
            
            if isempty(e)
                movefile(fullfile(p, n), fullfile(p, [n '.dcm']))
            end      
            
        end
              
    end
     
end

end