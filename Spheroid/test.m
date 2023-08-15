clear all;
close all
clc

sourcepath = 'C:\Users\Delmon\Desktop\Work\T47D';
destpath = 'C:\Users\Delmon\Dropbox\Jobb\Matlab\Spheroid';
folderList = getAllFolders(sourcepath);

struct = {};
sensitivity = [0.57 0.57 0.59 0.56 0.57 0.57 0.56];

for i = 1 % :length(folderList)
    
    [~, foldername, ~] = fileparts(folderList{i});
    fileList = getAllFiles(folderList{i});
    [~, ind] = sort(fileList);
    fileList = fileList(ind);
    
    filename_vec = {};
    A_outer_vec = [];
    radii_outer_vec = [];
    A_inner_vec = [];
    radii_inner_vec = [];
    
    counter = 0;
    
    for j = 1 % length(fileList) % 
        
        counter = counter + 1;
        
        [filepath, filename, ext] = fileparts(fileList{j});
        [img_RGB, img_MAP] = imread(fullfile(filepath, [filename ext]));
        
        warning('off', 'all')
        
        [A_outer, r_outer, A_inner, r_inner] = estimateSpheroidArea( ...
            img_RGB, img_MAP, foldername, filename, sensitivity(i));
        
%         struct.(sprintf('Days%s', foldername(1:2))). ...
%             (strrep(sprintf('%s', filename), '-', '_')).area = area;
%         struct.(sprintf('Days%s', foldername(1:2))). ...
%             (strrep(sprintf('%s', filename), '-', '_')).center = center;
%         struct.(sprintf('Days%s', foldername(1:2))). ...
%             (strrep(sprintf('%s', filename), '-', '_')).radii = radii;

        filename_vec{counter} = filename;
        A_outer_vec(counter) = A_outer;
        radii_outer_vec(counter) = r_outer;
        A_inner_vec(counter) = A_inner;
        radii_inner_vec(counter) = r_inner;
        
    end
        
%     header = {'Filename', 'Outer area (mm^2)', 'Outer radius (um)', ...
%         'Inner area (mm^2)', 'Inner radius (um)'};
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), header, 1, 'A1')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), filename_vec.', 1, 'A2')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), A_outer_vec.', 1, 'B2')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), radii_outer_vec.', 1, 'C2')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), A_inner_vec.', 1, 'D2')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), radii_inner_vec.', 1, 'E2')
    
    %close all
    
end
