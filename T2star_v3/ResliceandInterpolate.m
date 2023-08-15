function [] = ResliceandInterpolate(sourcepath, imgseq, ROIpatients)

folderList = getAllFolders(sourcepath);

%% Loop over each ROI for each patient
for i = 1:length(ROIpatients)
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        ROIpatients{i}.Patientname '\' imgseq]));
    
    if ~isempty(ind)
        
        ROIpatients{i}.Patientname
        folderList{ind}
        imgs = readDICOMMR(folderList{ind});
        slicelocationvec = [];
        temp_array = [];
        
        % Organize all 2D data matrices for every TE in a 4D array
        for j = 1:imgs.N_slices
            slicelocationvec = [slicelocationvec; imgs.slice{j}.slicelocation(1)];
            for k = 1:length(imgs.slice{1}.TE)
                temp_array(:,:,j,k) = imgs.slice{j}.data(:,:,k);
            end
        end
        
        % Compute new slice locations
        step = 5; % slice thickness (4 mm) + slice gap (1 mm) i DWI
        lower = (ROIpatients{i}.slicelocation - step):(step):slicelocationvec(end);
        upper = ROIpatients{i}.slicelocation:(-step):slicelocationvec(1);
        new_slicelocationvec = [lower(end:-1:1), upper].';
        
%         ROIpatients{i}.Patientname 
%         ROIpatients{i}.slicelocation
%         slicelocationvec
%         new_slicelocationvec
        
        % Reslice and interpolate
% %         [X,Y,Z] = meshgrid(1:size(temp_array,2), 1:size(temp_array,1), slicelocationvec);
% %         [X2,Y2,Z2] = meshgrid( ...
% %             1:ROIpatients{i}.PixelSpacing(2)/imgs.PixelSpacing(2):size(temp_array,2), ...
% %             1:ROIpatients{i}.PixelSpacing(1)/imgs.PixelSpacing(1):size(temp_array,1), ...
% %             new_slicelocationvec);
% %         for k = 1:size(temp_array,4)
% %             temp = temp_array(:,:,:,k);
% %             interpdata(:,:,:,k) = interp3(X, Y, Z, temp, X2, Y2, Z2, 'spline');
% %         end
        
        [X,Y] = meshgrid(1:size(EBT3_img,2), 1:size(EBT3_img,1));
        [X2,Y2] = meshgrid( ...
            1:Colony_img.PixelSpacing(2)/EBT3_img.PixelSpacing(2):size(EBT3_img,2), ...
            1:Colony_img.PixelSpacing(1)/EBT3_img.PixelSpacing(1):size(EBT3_img,1));
        interpdata(:,:) = interp2(X, Y, EBT3_img, X2, Y2, 'spline');
        
        % Write DICOM
        path = fullfile([sourcepath '\' ROIpatients{i}.Patientname], ...
            [imgseq 'reslice_v2']);
        if ~exist(path, 'dir')
            mkdir(path);
        end
        
        temp_info = imgs.slice{1}.info;
        
        for j = 1:size(interpdata, 3)
            
            temp_info.SeriesNumber = temp_info.SeriesNumber + 80;
            temp_info.SeriesDescription = [imgseq 'Reslice'];
            temp_info.PixelSpacing = [ROIpatients{i}.PixelSpacing(1) ...
                ROIpatients{i}.PixelSpacing(2)];
            temp_info.SliceLocation = new_slicelocationvec(j);
            temp_info.ImagePositionPatient(3) = new_slicelocationvec(j);
            
            for k = 1:size(interpdata, 4)
                temp_info.EchoTime = imgs.slice{1}.TE(k);
                temp_info.EchoNumbers = k;
                dicomwrite(uint16(interpdata(:,:,j,k)), ...
                    fullfile(path, strcat([imgseq 'Reslice'], ...
                    num2str(j + length(new_slicelocationvec)*(k - 1)), '.dcm')), ...
                    temp_info, 'CreateMode', 'Copy', 'CompressionMode', 'None');
%                 dlmwrite(fullfile(path, strcat([imgseq 'Reslice'], ...
%                     num2str(j + length(new_slicelocationvec)*(k - 1)), '.txt')), ...
%                     interpdata(:,:,j,k), 'precision', 16); 
            end
            
        end
        
                %         for k = 1:size(temp_array,4)
        %             temp = temp_array(:,:,:,k);
        %             temp = reshape(temp, [], size(temp_array,3)).';
        %             interp_images = interp1(slicelocationvec, temp, ...
        %                 new_slicelocationvec, 'pchip').';
        %             interp_images = reshape(interp_images, ...
        %                 [size(temp_array,1), size(temp_array,2), ...
        %                 length(new_slicelocationvec)]);
        %             final_interp_images(:,:,:,k) = interp_images;
        %         end
        
        %         % Plot
        %         dim = [0.2 0.5 0.3 0.3];
        %         for k = 1:size(temp_array, 4)
        %             data = temp_array(:,:,:,k);
        %             figure()
        %             for j = 1:size(temp_array, 3)
        %                 imshow(data(:,:,j), [min(data(:)) max(data(:))])
        %                 str = {['TE=', num2str(imgs.slice{1}.TE(k))], ...
        %                     ['SliceLocation=', num2str(new_slicelocationvec(j))]};
        %                 annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
        %                     'Color', 'white')
        %                 drawnow
        %                 pause(0.9)
        %             end
        %         end
           
    end
    
    clear imgs slicelocationvec temp_array step lower upper interpdata ...
    new_slicelocationvec X Y Z X2 Y2 Z2 temp path temp_info
    
end

clear ind

end