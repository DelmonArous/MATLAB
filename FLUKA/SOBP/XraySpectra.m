clear all;
close all;
fclose('all');
clc;

%% Directory
sourcepath  = 'C:\Users\delmo\Desktop\Spectra\DAT';
destpath    = 'C:\Users\delmo\Desktop\Spectra\TXT';
filelist    = getAllFiles(sourcepath);

%% Loop through X-ray spectra files
for i = 1:length(filelist)
    
    [p, fn, ~] = fileparts(filelist{i});
    
    % Read .dat file
    fileID = fopen(filelist{i}, 'r');
    
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    
    A = fscanf(fileID, formatSpec, sizeA);
    A = A.';
    
    fclose(fileID);
    
    e       = A(:,1)' ./ 10^6; % convert from kV to GV
    e_min   = [e(1)-mean(diff(e)) e(1:end-1)];
    e_max   = e;
    w       = A(:,2)';
    
    % Write .dat file
    fid = fopen(fullfile(destpath, [fn '.txt']), 'wt');
    fprintf(fid, '%.7f  %.7f   %.7f\n', [e_min; e_max; w]);
    fclose(fid);
    
end
