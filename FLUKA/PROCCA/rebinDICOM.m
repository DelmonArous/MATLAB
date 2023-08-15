clc
clear all
close all

sourcepath = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Mouse_male_20102020\Rec_DICOM 16 bits\DICOM';
destpath = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Mouse_male_20102020\Rec_DICOM 16 bits\DICOMreslice';
filelist = getAllFiles(sourcepath);

counter     = 0;
slice_step  = 12;
for i = 1:slice_step:length(filelist)
    
    [path, filename, ext] = fileparts(filelist{i});
    
    try
        % If dicominfo is successful, store the header information
        info = dicominfo(fullfile(path, filename), 'UseDictionaryVR', true);
        [~, MSGID] = lastwarn();
        warning('off', MSGID)
    catch
        exception
        if any(exception.identifier)
            try
                info = dicominfo(fullfile(path, filename), ...
                    'UseDictionaryVR', false);
                [~, MSGID] = lastwarn();
                warning('off', MSGID)
            catch
                warning(['File ', filename, ...
                    ' is not a valid DICOM object. Directory: ', path]);
                continue;
            end
        else
            % Otherwise, the file is either corrupt or not a real DICOM
            % file, so throw an error
            warning(['File ', filename, ...
                ' is not a valid DICOM object. Directory: ', path]);
            continue;
        end
    end
    
    img = double(dicomread(info));
    
    % Interpolate
    [X,Y] = meshgrid(1:size(img,2), 1:size(img,1));
%     [X2,Y2] = meshgrid(1:512, 1:512);
    [X2,Y2] = meshgrid(1:0.0503/info.PixelSpacing(2):size(img,2), ...
        1:0.05/info.PixelSpacing(1):size(img,1));
    interp_img = interp2(X, Y, img, X2, Y2, 'spline');
    
    % Adapt DICOM metadata
    temp_info                       = info;
%     temp_info.SeriesNumber          = temp_info.SeriesNumber + 80;
    temp_info.PatientID             = 'C57BL/6';
    temp_info.PatientName           = 'Mouse';
    temp_info.PatientSex            = 'M';
    temp_info.PatientAge            = '1';
    temp_info.InstitutionName       = 'Radiumhospital';
    temp_info.StudyDescription      = 'MicroCT_HeadAndNeck';
    temp_info.SeriesDescription     = 'MicroCT';
    temp_info.PixelSpacing          = [0.05 0.0503];
    temp_info.SliceThickness        = info.SliceThickness * slice_step;
    temp_info.InstanceNumber        = counter + 1;
    temp_info.FrameOfReferenceUID   = '1.3.6.1.4.1.14519.5.2.1.5099.8010.812355329601144442952825387096';
    temp_info.PositionReferenceIndicator = '';
    temp_info.StorageMediaFileSetUID = '1.3.6.1.4.1.14519.5.2.1.5099.8010.784780811979408194740251235264';
    temp_info.WindowCenter          = [35; 700];
    temp_info.WindowCenterWidthExplanation = 'WINDOW1\WINDOW2';
    temp_info.WindowWidth           = [80; 3200];
    
%     figure(); imshow(interp_img, [])
%     [info.RescaleSlope info.RescaleIntercept]
%     [min(interp_img(:)) max(interp_img(:))]
    interp_img = ((interp_img - 14000) ./ (14000 - 6627)) .* 1000;
    interp_img = interp_img + 1024;
    interp_img(interp_img <= 0) = 0; 
%     interp_img.*info.RescaleSlope + info.RescaleIntercept;
%     [min(interp_img(:)) max(interp_img(:))]
%     [min(int16(interp_img(:))) max(int16(interp_img(:)))]
%     [min(uint16(interp_img(:))) max(uint16(interp_img(:)))]
%     figure(); imshow(interp_img, [])

    % Write DICOM file
    dicomwrite(uint16(interp_img), ...
        fullfile(destpath, [sprintf('_%04d', counter) ext]), ...
        temp_info, 'CreateMode', 'Copy', 'CompressionMode', 'None', ...
        'Endian', 'ieee-le');
    
    %     % Copy file
    %     copyfile(filelist{i}, fullfile(destpath, ...
    %         [sprintf('%04d', counter) ext]))
    
    counter = counter + 1;
    
end

img = readDICOMimage(...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Mouse_male_20102020\Rec_DICOM 16 bits\DICOMreslice');
noe = writeDICOMimage(img, ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\temp\Mouse_male_20102020\Rec_DICOM 16 bits\DICOMreslice_v2');


