clear all;
close all
clc

sourcepath = 'C:\Users\Delmon\Desktop\Work\T47D';
destpath = 'C:\Users\Delmon\Dropbox\Jobb\Matlab\Spheroid';
folderList = getAllFolders(sourcepath);

struct = {};

for i = 7 % :length(folderList)
    
    [~, foldername, ~] = fileparts(folderList{i});
    fileList = getAllFiles(folderList{i});
    [~, ind] = sort(fileList);
    fileList = fileList(ind);
        
    filename_vec = {};
    radii_vec = [];
    area_vec = [];
    
    counter = 0;
 
    for j = 1:28 % 19:28 % 1:12 % 13:16 % 1:length(fileList)
        
        counter = counter + 1;
        
        [filepath, filename, ext] = fileparts(fileList{j});
        [img_RGB, img_MAP] = imread(fullfile(filepath, [filename ext]));
        
        warning('off', 'all')
        
        [area, center, radii] = estimateSpheroidNecroticArea( ...
            img_RGB, img_MAP, foldername, filename);
        
%         struct.(sprintf('Days%s', foldername(1:2))). ...
%             (strrep(sprintf('%s', filename), '-', '_')).area = area;
%         struct.(sprintf('Days%s', foldername(1:2))). ...
%             (strrep(sprintf('%s', filename), '-', '_')).center = center;
%         struct.(sprintf('Days%s', foldername(1:2))). ...
%             (strrep(sprintf('%s', filename), '-', '_')).radii = radii;
        
%         filename_vec{counter} = filename;
%         radii_vec(counter) = radii;
%         area_vec(counter) = area;
%         
    end
        
%     header = {'Filename', 'Radius (um)', 'Area (mm^2)'};
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), header, 1, 'A1')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), filename_vec.', 1, 'A2')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), radii_vec.', 1, 'B2')
%     xlswrite(fullfile(destpath, [foldername '.xlsx']), area_vec.', 1, 'C2')

    % close all
    
end
