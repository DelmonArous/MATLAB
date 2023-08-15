function [ROIpatients] = computeR2starROI(sourcepath, ROIpatients)

folderList = getAllFolders(sourcepath);

%% Cohort 1 + Cohort 2
xvec = [1 1 24 0 20 5 19 -8 -1 4 -8 17 -3 23 -3 ...
    -5 10 9 0 4 4 22 -6 -19 ...
    23 -1 10 13 24 8 14 15 11 ...
    26 13 -5 9 6 10 10 -6 13 ...
    4 2 25 11 ...
    1 14 19 2 12 8 4 8 ... %% Cohort 2
    -9 9 12 29 10 16 17 3 27 ...
    1 2 11 25 1 2 5 2 4 ...
    2 24 3 4 6 12 -6 28 12 ...
    9 10 10 6 7 6 2 4 33 4 ...
    13 -11 -7 -11 11 20];
yvec = [1 1 -23 1 -18 -4 -18 9 2 -3 10 -17 5 -21 5 ...
    6 -10 -7 2 -3 -2 -20 7 21 ...
    -22 2 -8 -10 -22 -6 -14 -13 -10 ...
    -24 -12 6 -7 -6 -8 -9 7 -12 ...
    -3 -1 -25 -10 ...
    1 -13 -18 -2 -10 -7 -2 -6 ... %% Cohort 2
    9 -9 -10 -28 -8 -16 -16 -2 -27 ...
    0 -1 -10 -24 -1 -1 -4 -2 -2 ...
    0 -23 -2 -4 -4 -11 7 -26 -11 ...
    -8 -10 -9 -4 -6 -4 -2 -2 -33 -2 ...
    -11 12 8 12 -11 -18];

%% Loop over each ROI for each patient
for i = 1:length(ROIpatients)
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        ROIpatients{i}.Patientname '\R2starreslice_v2']));
    
    if ~isempty(ind)
        
        fileList = getAllFiles(folderList{ind});
        
        for j = 1:numel(fileList)
            
            [p, n, e] = fileparts(fileList{j});
            
            %% Read only T2* maps saved as DICOM
            if (contains(n, 'T2starMap') || contains(n, 'R2starReslice')) ...
                    && strcmp(e, '.dcm')
                
                try
                    % If dicominfo is successful, store the header information
                    info = dicominfo(fileList{j}, 'UseDictionaryVR', true);
                catch
                    % Otherwise, the file is either corrupt or not a real DICOM
                    % file, so throw an error
                    warning(['File ', n, ' is not a valid DICOM object.' ...
                        newline 'Directory: ', p]);
                    continue;
                end
                  
                % Check if stored DICOM slice corresponds to the ROI slice -
                % do computation only on the stored DICOM T2* slice
                if abs((info.ImagePositionPatient(3) - ROIpatients{i}.slicelocation) / ...
                        info.ImagePositionPatient(3)) <= 0.0001
                    
%                      ROIpatients{i}.Patientname
%                     [info.ImagePositionPatient(3) ROIpatients{i}.slicelocation]
                    
                    T2starMap = double(dicomread(info)) ./ 100;
                    R2starMap = (1./T2starMap) .* 1000; % in s^(-1)
                    
                    R2starpxX = abs(ROIpatients{i}.mmX - ...
                        info.ImagePositionPatient(1)) ./ info.PixelSpacing(1);
                    R2starpxY = abs(ROIpatients{i}.mmY - ...
                        info.ImagePositionPatient(2)) ./ info.PixelSpacing(2);
                    R2starROImask = double(poly2mask(R2starpxX, R2starpxY, ...
                        double(info.Rows), double(info.Columns)));
                    
                    %% Store information
                    ROIpatients{i}.R2starMap = R2starMap;
                    ROIpatients{i}.R2starROImask = R2starROImask;
                    ROIpatients{i}.R2starImagePositionPatient = info.ImagePositionPatient;
                    
                    %%
                    tempmask = zeros(size(ROIpatients{i}.R2starMap));
                    s = size(ROIpatients{i}.ADCROImask);

                    x = (ROIpatients{i}.ImagePositionPatient(1) - ...
                        info.ImagePositionPatient(1)) / ...
                        ROIpatients{i}.PixelSpacing(1);
                    y = (ROIpatients{i}.ImagePositionPatient(2) - ...
                        info.ImagePositionPatient(2)) / ...
                        ROIpatients{i}.PixelSpacing(2);
                    
                    x = x + xvec(i); y = y + yvec(i);
                    tempmask(x:x+s(1)-1, y:y+s(2)-1) = ROIpatients{i}.ADCROImask;
                    ROIpatients{i}.R2starROImask = tempmask;
 
%                     h = figure();
%                     imshow(translatedADCROImask, [0 1])
%                     title([ROIpatients{i}.Patientname ' ADC ROI translated'])
% 
%                     h = figure();
%                     imshow(ROIpatients{i}.R2starROImask, [0 1])
%                     title([ROIpatients{i}.Patientname ' R2* ROI'])

                    %%
                    % Multiply the R2* array by the binary structure mask and reshape
                    % into a 1D vector (to keep zero values within the
                    % mask of the given structure, 1e-6 is added)
                    ROIpatients{i}.ROIR2starMap = ROIpatients{i}.R2starMap .* ...
                        ROIpatients{i}.R2starROImask;
                    ROIR2starvalues = reshape((ROIpatients{i}.R2starMap + 1e-6) .* ...
                        ROIpatients{i}.R2starROImask, 1, []);
                    
                    % Remove voxel values inside the mask of the structure,
                    % in which are zero or NaN
                    ROIR2starvalues(ROIR2starvalues == 0) = [];
                    ROIR2starvalues(isnan(ROIR2starvalues)) = [];
                    ROIpatients{i}.ROIR2starvalues = ROIR2starvalues;
                    
                    % Calculate the 1-99th R2* percentiles
                    ROIpatients{i}.R2starmedian = median(ROIR2starvalues);
                    ROIpatients{i}.R2starpercentile = prctile(ROIR2starvalues, 1:99);
                                       
                    %% Plot
%                     h = figure();
%                     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%                     imshow(ROIpatients{i}.R2starMap, [0 50], 'colormap', jet(256))
%                     % [min(ROIpatients{i}.R2starMap(:)) max(ROIpatients{i}.R2starMap(:))], ...
%                     c = colorbar;
%                     c.Label.String = 'R2* (s^{-1})';
%                     hold on
%                     contour(ROIpatients{i}.R2starROImask, 1, 'Linestyle', '-.', ...
%                         'LineColor', 'black', 'LineWidth', 2.0)
%                     leg = legend('Tumor outline');
%                     title(leg, ROIpatients{i}.Patientname)
%                     set(gca, 'FontSize', 16)
%                     shading interp
%                     hold off
%                     saveas(h, sprintf('R2starMapWithROI_%s.png', ...
%                         ROIpatients{i}.Patientname));
                    
                end
                
            end
            
        end
        
    end
    
end

clear folderList fileList ind p n e info R2starMap R2starROImask ROIvalues

end