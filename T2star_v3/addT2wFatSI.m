function [R2starpatients] = addT2wFatSI(T2fatpath, R2starpatients) 

%% Get all T2w fat SI files
fileList_fatSI = getAllFiles(T2fatpath);

%%
for i = 1:length(R2starpatients)
    
    ind = find(strcmp(fileList_fatSI, ...
        horzcat(T2fatpath, '\', R2starpatients{i}.Patientname, '.txt')));
    
    if ~isempty(ind)
       % Read intensity of T2w fat signal
        R2starpatients{i}.SIfat = double(dlmread(fileList_fatSI{ind})); 
    end
    
end

end