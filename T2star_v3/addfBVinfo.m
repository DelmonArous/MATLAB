function [ROIpatients] = addfBVinfo(sourcepath, ROIpatients)

folderList = getAllFolders(sourcepath);

for i = 1:length(ROIpatients)
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        ROIpatients{i}.Patientname '\fBV']));
    
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
                
                %% Compute fBV map
                ROIpatients{i}.fBVmap = double(dicomread(info)) ./ 10^4;
                
                % Multiply the fBV array by the binary structure mask and reshape
                % into a 1D vector (to keep zero values within the
                % mask of the given structure, 1e-6 is added)
                ROIpatients{i}.ROIfBVmap = ROIpatients{i}.fBVmap .* ...
                    ROIpatients{i}.ADCROImask;
                ROIfBVvalues = reshape((ROIpatients{i}.fBVmap + 1e-6) .* ...
                    ROIpatients{i}.ADCROImask, 1, []);
                
                % Remove voxel values inside the mask of the structure,
                % in which are zero or NaN
                ROIfBVvalues(ROIfBVvalues == 0) = [];
                ROIfBVvalues(isnan(ROIfBVvalues)) = [];
                ROIpatients{i}.ROIfBVvalues = ROIfBVvalues;
                
                % Calculate the 1-99th fBV percentiles
                ROIpatients{i}.fBVmedian = median(ROIfBVvalues);
                ROIpatients{i}.fBVpercentile = prctile(ROIfBVvalues, 1:99);
                
                %% Plot              
%                 h = figure();
%                 set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%                 imshow(ROIpatients{i}.fBVmap, ...
%                     [min(ROIpatients{i}.fBVmap(:)) 0.9], ... % max(ROIpatients{i}.fBVmap(:))], ...
%                     'colormap', jet(256))
%                 c = colorbar;
%                 c.Label.String = 'fBV (a.u)';
%                 set(gca, 'FontSize', 16)
%                 hold on
%                 contour(ROIpatients{i}.ADCROImask, 1, 'Linestyle', '-.', ...
%                     'LineColor', 'black', 'LineWidth', 2.0)
%                 leg = legend('Tumor outline');
%                 title(leg, ROIpatients{i}.Patientname)
%                 shading interp
%                 hold off
%                 saveas(h, sprintf('fBVMapWithROI_%s.png', ...
%                     ROIpatients{i}.Patientname));
                
            end
            
            clear info ROIfBVvalues
            
        end
        
        clear p n
        
    end
    
    clear fileList
    
end

end