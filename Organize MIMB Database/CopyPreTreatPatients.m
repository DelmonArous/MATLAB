function [] = CopyPreTreatPatients(sourcepath, destpath, xlsxdata)

Nr = cell2mat(xlsxdata(2:end-1, 1));
MedInsightID = cell2mat(xlsxdata(2:end-1, 2));
TreatStart = cell2mat(xlsxdata(2:end-1, 3));
TreatStart = TreatStart(:, 1:10);
patIDvec = strings(length(Nr), 1);

for j = 1:length(Nr)
    patIDvec(j) = strcat('Pat', num2str(Nr(j)), '_', num2str(MedInsightID(j)));
end

folderList = getAllFolders(sourcepath);
for i = 1:length(folderList)
    
    [folderpath, foldername, ~] = fileparts(folderList{i});
    
    for j = 1:length(patIDvec)
    
        if strcmp(foldername, patIDvec(j))
            foldername
            initiateCopyPreTreatPatient(fullfile(folderpath, foldername), ...
                destpath, TreatStart(j,1:10));
        end
    
    end
end

end