function [ROIpatients] = addADCinfo(sourcepath, ROIpatients)

folderList = getAllFolders(sourcepath);

for i = 1:length(ROIpatients)
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        ROIpatients{i}.Patientname '\ADC_v2']));
    
    if ~isempty(ind)
        
        fileList = getAllFiles(folderList{ind});
        
        for j = 1:length(fileList)
            
            [p, n, ~] = fileparts(fileList{j});
            
            if strcmp(ROIpatients{i}.Filename(end-1:end), n(end-1:end))
                
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
                
                ROIpatients{i}.Patientname
                
                %% Store ADC slice metadata information
                ROIpatients{i}.info = info;
                ROIpatients{i}.slicelocation = info.SliceLocation;
                ROIpatients{i}.PixelSpacing = info.PixelSpacing;
                ROIpatients{i}.ImagePositionPatient = info.ImagePositionPatient;
                ROIpatients{i}.mmX = (ROIpatients{i}.pxX .* ...
                    info.PixelSpacing(1)) + info.ImagePositionPatient(1);
                ROIpatients{i}.mmY = (ROIpatients{i}.pxY .* ...
                    info.PixelSpacing(2)) + info.ImagePositionPatient(2);
                
                %% Compute ADC map and ADC-ROI mask, and store the polygon
                ADCROImask = double(poly2mask(ROIpatients{i}.pxX, ...
                    ROIpatients{i}.pxY, double(info.Rows), double(info.Columns)));
                ROIpatients{i}.ADCROImask = ADCROImask;
                ROIpatients{i}.ADCmap = double(dicomread(info))./ 10^6;
                
                % Multiply the ADC array by the binary structure mask and reshape
                % into a 1D vector (to keep zero values within the
                % mask of the given structure, 1e-6 is added)
                ROIpatients{i}.ROIADCmap = ROIpatients{i}.ADCmap .* ...
                    ROIpatients{i}.ADCROImask;
                ROIADCvalues = reshape((ROIpatients{i}.ADCmap + 1e-6) .* ...
                    ROIpatients{i}.ADCROImask, 1, []);
                
                % Remove voxel values in which are zero or NaN
                %(basically, voxels outside of the structure mask
                ROIADCvalues(ROIADCvalues == 0) = [];
                ROIADCvalues(isnan(ROIADCvalues)) = [];
                ROIpatients{i}.ROIADCvalues = ROIADCvalues;
                
                % Calculate the 1-99th ADC percentiles
                ROIpatients{i}.ADCmedian = median(ROIADCvalues);
                ROIpatients{i}.ADCpercentile = prctile(ROIADCvalues, 1:99);
                
                %% Plot      
%                 h = figure();
%                 set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%                 imshow(ROIpatients{i}.ADCmap, ...
%                     [min(ROIpatients{i}.ADCmap(:)) max(ROIpatients{i}.ADCmap(:))], ...
%                     'colormap', jet(256))
%                 c = colorbar;
%                 c.Label.String = 'ADC (mm^2/s)';
%                 set(gca, 'FontSize', 16)
%                 hold on
%                 contour(ADCROImask, 1, 'Linestyle', '-.', ...
%                     'LineColor', 'black', 'LineWidth', 2.0)
%                 leg = legend('Tumor outline');
%                 title(leg, ROIpatients{i}.Patientname)
%                 shading interp
%                 hold off
%                 saveas(h, sprintf('ADCMapWithROI_%s.png', ...
%                     ROIpatients{i}.Patientname));
                
            end
            
            clear info ADCROImask ROIADCvalues
            
        end
        
        clear p n
        
    end
    
    clear fileList
    
end

end