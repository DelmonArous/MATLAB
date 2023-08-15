
% << How many [%] of the image the given color occupies >>
%
% Input: you have to know the RGB bands of the color in which you are interested
%       The input RGB can be inserted as bands (so that slight
%       fluctuations of the colors are also captured) or exact number
% Output: Several images are generated during the analysis:
%           Red, Blue and Green masks of the original color image
%           Histograms of red, blue, and green parts of the image
%           Distribution of the size of the blobs
%           Mask excluding the small blobs
%           Mask with the filled holes
%           Comaprison of the original and final image
%           !!! Images are not saved !!!
%         Table of the summary - number of the blobs, their areas [pixel]
%           and their color
%           The resultant table is saved in xls file - blobs' area, color,
%               color bands (user input), minimal size of the blobs (user
%               input), % of the image coverd by the chosen color
% Notes:
%  1) you can desides wheather you want to get rid of small blobls of
%       the color areas (area smaller than given "user-input" value
%       will be excluded)
%  2) you are asked if you want to fill the holes in the blobs found.
%  3) use imtool to explore the RBG colors in the image
% ------------------------------
% 24.2.2014
% Michala.Cadova@zzm.uzh.ch
% Laboratory of Biomechanics, Czech Technical University in Prague, Czech Republic
% Michala.Cadova@fs.cvut.cz
% Dental Clinic, University of Zurich, Switzerland
% Michala.Cadova@zzm.uzh.ch
% ------------------------------
% Literature:
%  1) rgb table: http://www.rapidtables.com/web/color/RGB_Color.htm
%  2) matlab webinar: http://www.mathworks.ch/videos/medical-imaging-workflows-with-matlab-81850.html?form_seq=conf924&confirmation_page&wfsid=5335062
%  3) SimpleColorDetection.m - http://www.mathworks.ch/matlabcentral/fileexchange/26420-simplecolordetection/all_files
%  4) bwareaopen - http://www.mathworks.ch/ch/help/images/ref/bwareaopen.html
% ------------------------------
%         CHANGES
% ------------------------------
% 27.2.2014 - V13
%  - Create new folder for the results (name based on the image name)
%  - Automatically save images
%  - xls file with the resutls saved in the result folder
%  - headers in the xls file
%  - saving the xls file - decission window
%  - histogram of the distribution of the size of the blobs
%  - new image of the result
% ------------------------------
clear all; close all; clc
%% General settings
fontSize = 14;
%% Pick the image and Load the image
[filename, pathname] = uigetfile( ...
    {'*.jpg;*.tif;*.tiff;*.png;*.bmp', 'All image Files (*.jpg, *.tif, *.tiff, *.png, *.bmp)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Pick a file');
f = fullfile(pathname, filename);
disp('Reading image')
rgbImage = imread(f);
[rows columns numberOfColorBands] = size(rgbImage);
%% Create directory for the results
k = findstr('.', filename); % find where the extension starts
FileToSave =  filename(1:k-1); %
fileNameWrongCounter  = 0;
% if the name of the input wav file is longer than 256 chrachters,
% which is the limit for windows folder name's length,
% the directory will be composed of the first 200 chraracters
% of the file + the "CounterNumber" will be added to see,
% how meny files has to long name
% However, because of the excel file names restriction, the file name
% cannot be 256, but only 218. Therefore, the folder name lengthe is
% restricted to 150
CreateNewFoldrSucces = 0;
while CreateNewFoldrSucces == 0
    if length(filename) <  150
        % the lenght of the name for the folder can be up to 255
        % characters, but later, the lenght of the file name (TOGETHER WITH
        % the path) of the excel/csv files with the results can be only 218
        DirectoryName = ['Res_', FileToSave];
        % Matalb does not like when the names beginns with a NUMBER!!!
        % F = folder
        % dont apend the extension of the wav file: (.wav) = 4characters
        % FileToSave is already adjusted with respect to the
        % prescribed file name, and is without the extension
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(DirectoryName);
        % SUCCESS = 1 -> folder was created.
        % what if the folder already exists?
    else % if the input file name is longer than 150 characters
        fileNameWrongCounter = fileNameWrongCounter + 1;
        DirectoryName = strcat('Res_', FileToSave(1:30), '_', num2str(fileNameWrongCounter));
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(DirectoryName );
    end
    
    if SUCCESS == 1
        CreateNewFoldrSucces = 1;
        break
    end
end

%% If the image is monochrome (indexed), convert it to color.
% Check to see if it's an 8-bit image needed later for scaling).
if strcmpi(class(rgbImage), 'uint8')
    % Flag for 256 gray levels.
    eightBit = true;
else
    eightBit = false;
end

if numberOfColorBands == 1
    if isempty(storedColorMap)
        % Just a simple gray level image, not indexed with a stored color map.
        % Create a 3D true color image where we copy the monochrome image into all 3 (R, G, & B) color planes.
        rgbImage = cat(3, rgbImage, rgbImage, rgbImage);
    else
        % It's an indexed image.
        rgbImage = ind2rgb(rgbImage, storedColorMap);
        % ind2rgb() will convert it to double and normalize it to the range 0-1.
        % Convert back to uint8 in the range 0-255, if needed.
        if eightBit
            rgbImage = uint8(255 * rgbImage);
        end
    end
end

%% Display the color image
disp('Displaying color original image')
F1 = figure(1);
subplot(3,4,1);
imshow(rgbImage);

if numberOfColorBands > 1
    title('Original Color Image', 'FontSize', fontSize);
else
    caption = sprintf('Original Indexed Image\n(converted to true color with its stored colormap)');
    title(caption, 'FontSize', fontSize);
end

%% Size of the picture - to occupy the whole screen
scnsize = get(0,'ScreenSize'); % - - width height
position = get(F1,'Position'); % x-pos y-pos widht height
outerpos = get(F1,'OuterPosition');
borders = outerpos - position;
edge = abs(borders(1))/2;
pos1 = [edge,...
    1/20*scnsize(4), ...
    9/10*scnsize(3),...
    9/10*scnsize(4)];
set(F1,'OuterPosition',pos1)
%% Explore RGB
% Extract out the color bands from the original image
% into 3 separate 2D arrays, one for each color component.
redBand = rgbImage(:, :, 1);
greenBand = rgbImage(:, :, 2);
blueBand = rgbImage(:, :, 3);
% Display them.
subplot(3, 4, 2);
imshow(redBand);
title('Red Band', 'FontSize', fontSize);
subplot(3, 4, 3);
imshow(greenBand);
title('Green Band', 'FontSize', fontSize);
subplot(3, 4, 4);
imshow(blueBand);
title('Blue Band', 'FontSize', fontSize);
%% Compute and plot the red histogram.
hR = subplot(3, 4, 6);
[countsR, grayLevelsR] = imhist(redBand);
maxGLValueR = find(countsR > 0, 1, 'last');
maxCountR = max(countsR);
bar(countsR, 'r');
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
title('Histogram of Red Band', 'FontSize', fontSize);
%% Compute and plot the green histogram.
hG = subplot(3, 4, 7);
[countsG, grayLevelsG] = imhist(greenBand);
maxGLValueG = find(countsG > 0, 1, 'last');
maxCountG = max(countsG);
bar(countsG, 'g', 'BarWidth', 0.95);
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
title('Histogram of Green Band', 'FontSize', fontSize);
%% Compute and plot the blue histogram.
hB = subplot(3, 4, 8);
[countsB, grayLevelsB] = imhist(blueBand);
maxGLValueB = find(countsB > 0, 1, 'last');
maxCountB = max(countsB);
bar(countsB, 'b');
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
title('Histogram of Blue Band', 'FontSize', fontSize);
%% Set all axes to be the same width and height.
% This makes it easier to compare them.
maxGL = max([maxGLValueR,  maxGLValueG, maxGLValueB]);
if eightBit
    maxGL = 255;
end
maxCount = max([maxCountR,  maxCountG, maxCountB]);
axis([hR hG hB], [0 maxGL 0 maxCount]);
%% Plot all 3 histograms in one plot.
subplot(3, 4, 5);
plot(grayLevelsR, countsR, 'r', 'LineWidth', 2);
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
hold on;
plot(grayLevelsG, countsG, 'g', 'LineWidth', 2);
plot(grayLevelsB, countsB, 'b', 'LineWidth', 2);
title('Histogram of All Bands', 'FontSize', fontSize);
maxGrayLevel = max([maxGLValueR, maxGLValueG, maxGLValueB]);
% Trim x-axis to just the max gray level on the bright end.
if eightBit
    xlim([0 255]);
else
    xlim([0 maxGrayLevel]);
end
%% Now select thresholds for the 3 color bands.
% pop-up window
prompt = {'RED color threshold LOWER:','RED color threshold UPPER:', ...
    'GREEN color threshold LOWER:','GREEN color threshold UPPER:', ...
    'BLUE color threshold LOWER:','BLUE color threshold UPPER:'};
dlg_title = 'Input';
num_lines = 1;
% def = {'0','50', '150','255','0','100'}; % only green
def = {'45','150', '50','150','0','50'}; % green and similar
% def = {'235','255','235','255','0','200'};% red
% def = {'250','255', '250','255','0','245'};% yellow
% def = {'190','255', '160','245','60','165'};% yellow + light brown + yellow-white
answer = inputdlg(prompt,dlg_title,num_lines,def);

redThresholdLow = str2num(answer{1});
redThresholdHigh = str2num(answer{2});
greenThresholdLow = str2num(answer{3});
greenThresholdHigh = str2num(answer{4});
blueThresholdLow = str2num(answer{5});
blueThresholdHigh = str2num(answer{6});

%% Show the thresholds as vertical red bars on the histograms.
PlaceThresholdBars(1, 3,4, 6, redThresholdLow, redThresholdHigh, fontSize, max(countsR));
PlaceThresholdBars(1, 3,4, 7, greenThresholdLow, greenThresholdHigh,fontSize, max(countsG));
PlaceThresholdBars(1, 3,4, 8, blueThresholdLow, blueThresholdHigh,fontSize, max(countsB));

%% Now apply each color band's particular thresholds to the color band
redMask = (redBand >= redThresholdLow) & (redBand <= redThresholdHigh);
greenMask = (greenBand >= greenThresholdLow) & (greenBand <= greenThresholdHigh);
blueMask = (blueBand >= blueThresholdLow) & (blueBand <= blueThresholdHigh);

%% Display the thresholded binary images.
subplot(3, 4, 10);
imshow(redMask, []);
title('Is-Red Mask', 'FontSize', fontSize);
subplot(3, 4, 11);
imshow(greenMask, []);
title('Is-Green Mask', 'FontSize', fontSize);
subplot(3, 4, 12);
imshow(blueMask, []);
title('Is-Blue Mask', 'FontSize', fontSize);

%% Combine the masks to find where all 3 are "true."
% Then we will have the mask of only the chosen color parts of the image.
ObjectsMask = uint8(redMask & greenMask & blueMask);
subplot(3, 4, 9);
imshow(ObjectsMask, []);
caption = sprintf('Mask of the objects with chosen color');
title(caption, 'FontSize', fontSize);

f2 = fullfile(pathname, DirectoryName);
filNamePlot =  strcat(f2, '\Figure_1.png');
%saveas(gcf,num2str(char(filNamePlot)))
%% Histogram small areas
% Measure the mean RGB and area of all the detected blobs.
[meanRGB, areas, numberOfBlobs] = MeasureBlobs(ObjectsMask, redBand, greenBand, blueBand);
F30 = figure(30);
plot(areas(:,1))
title('Distribution of the areas of the blobls')

save('BlobAreas', 'meanRGB', 'areas', 'numberOfBlobs')
% Size of the picture - to occupy the whole screen
pos2 = [1/4*scnsize(3),...
    1/20*scnsize(4), ...
    2/3*scnsize(3),...
    2/3*scnsize(4)];
set(F30,'OuterPosition',pos2)
xlabel('Number of blobs/"islands"')
ylabel('Area of the blobs [pixels]')
filNamePlot =  strcat(f2, '\Figure_30.png');
%saveas(gcf,num2str(char(filNamePlot)))

figure(31)
XTickDescr = [0,5,10,20,50,100,200,300,500,1000,2000,3000];
N = hist(areas(:,1), XTickDescr);
bar(N)
set(gca,'XTick',1:length(XTickDescr))
set(gca,'XTickLabel',XTickDescr)
xlhand = get(gca,'xlabel');
set(xlhand,'string','X','fontsize',0.3)
xlabel('Size of the blobs [pixel] (each bar includes current value, exlude value to the right)')
ylabel('Number of blobs of current size')
zl = max(areas(:,1));
filNamePlot =  strcat(f2, '\Figure_31_histBlobs.png');
%saveas(gcf,num2str(char(filNamePlot)))
%% Insert the minimal area that will be counted
% Every blob smaller than this one will be ommited
prompt = {'Minimal area [pixels] of the blob that will be kept:'};
dlg_title = 'Min.area';
num_lines = 1;
def = {'10'};% yellow
answer2 = inputdlg(prompt,dlg_title,num_lines,def);
answer2 = str2num(answer2{1});

%% Ignore all small areas
F70 = figure(70);
subplot(2,2,1)
imshow(ObjectsMask, []);
title('Original mask', 'FontSize', fontSize)
set(F70,'OuterPosition',pos1)

ObjectsMask = uint8(bwareaopen(ObjectsMask,answer2));
figure(70)
subplot(2,2,2)
imshow(ObjectsMask, []);
title('Only the big particles', 'FontSize', fontSize)
%% Fill in any holes in the regions, since they are most likely red also.
message = sprintf('Do you want to close the holes in the blobs?');
reply = questdlg(message, 'Close holes?', 'Yes','No', 'Yes');
if strcmpi(reply, 'Yes')
    %figure(70)
    %     subplot(1,2,1);
    %         imshow(ObjectsMask, []);
    %         title('Original mask', 'FontSize', fontSize)
    subplot(2,2,3);
    ObjectsMask = uint8(imfill(ObjectsMask, 'holes'));
    imshow(ObjectsMask, []);
    title('Mask with filled holes', 'FontSize', fontSize);
end

filNamePlot =  strcat(f2, '\Figure_70.png');
%saveas(gcf,num2str(char(filNamePlot)))
%% Use the object mask to mask out the portions of the rgb image.
maskedImageR = ObjectsMask .* redBand;
maskedImageG = ObjectsMask .* greenBand;
maskedImageB = ObjectsMask .* blueBand;
% Concatenate the masked color bands to form the rgb image.
maskedRGBImage = cat(3, maskedImageR, maskedImageG, maskedImageB);

%% Show the masked off and original image.
F100 = figure(100);
subplot(1,2,1);
imshow(maskedRGBImage);
caption = sprintf('Masked Original Image');
title(caption, 'FontSize', fontSize);
subplot(1,2,2);
imshow(rgbImage);
title('The Original Image', 'FontSize', fontSize);

set(F100,'OuterPosition',pos1)

Text_descr =  strcat('R_l =', num2str(redThresholdLow), ', R_u = ', num2str(redThresholdHigh), ...
    ', G_l =', num2str(greenThresholdLow), ', G_u = ', num2str(greenThresholdHigh), ...
    ', B_l =', num2str(blueThresholdLow), ', B_u = ', num2str(blueThresholdHigh))
Text_descr2 =  strcat('Filled holes = ', reply, ', minimal counted area = ', num2str(answer2), ' [pixel]', ...
    ', max. blob size = ', num2str(max(areas(:,1))), ' [pixel], # of counted blobs = ', num2str(length(areas(:,1))));

text(-3000,2200,Text_descr) ;
text(-3000,2500,Text_descr2) ;

f2 = fullfile(pathname, DirectoryName);
filNamePlot =  strcat(f2, '\Figure_100.png');
%saveas(gcf,num2str(char(filNamePlot)))

%% Measure the mean RGB and area of all the detected blobs.
clear meanRGB
clear areas
clear numberOfBlobs
[meanRGB, areas, numberOfBlobs] = MeasureBlobs(ObjectsMask, redBand, greenBand, blueBand);
if numberOfBlobs > 0
    fprintf(1, '\n----------------------------------------------\n');
    fprintf(1, 'Blob #, Area in Pixels, Mean R, Mean G, Mean B\n');
    fprintf(1, '----------------------------------------------\n');
    for blobNumber = 1 : numberOfBlobs
        fprintf(1, '#%5d, %14d, %6.2f, %6.2f, %6.2f\n', blobNumber, areas(blobNumber), ...
            meanRGB(blobNumber, 1), meanRGB(blobNumber, 2), meanRGB(blobNumber, 3));
    end
else
    % Alert user that no  blobs were found.
    uiwait(msgbox('No blobs of given color were found in the image'), 'Error', 'error')
end
%% Compute how many % of the image the restult takes
OrigImageArea = rows*columns % pxls ... 100%
AreaOfChosenColor = sum(areas(:,1)) % ... x %
Area_Procent = (AreaOfChosenColor/OrigImageArea)*100
uiwait(msgbox(sprintf('%s %0.5g %s','The chosen color covers ', Area_Procent, ' % of the whole image.'), 'Results', 'modal'))
%% Plot only the resulted masked image - change the mask color to white
F100 = figure(101);
imshow(maskedRGBImage);
% caption = sprintf(['Masked Original Image. Given color covers ', num2str(Area_Procent), ' % of the area']);
% title(caption, 'FontSize', fontSize);
T101 = title(['Masked Original Image. Given color covers ', num2str(Area_Procent), ' % of the area.']);
set(T101, 'FontSize', fontSize);
set(F100,'OuterPosition',pos1)

f2 = fullfile(pathname, DirectoryName);
filNamePlot =  strcat(f2, '\Figure_101_Res.png');
%saveas(gcf,num2str(char(filNamePlot)))

%% save into xls file
message = sprintf('Do you want to save the resuts?');
reply = questdlg(message, 'Save results?', 'Yes','No', 'Yes');
if strcmpi(reply, 'Yes')
    
    XlsFilesave = fullfile(pathname, DirectoryName, 'BlobsResults');
    RGBHeadings = {'Red band','Green band' , 'Blue band'};
    
    xlswrite(XlsFilesave, RGBHeadings,'meanRGB', 'A1') ;
    xlswrite(XlsFilesave, meanRGB,'meanRGB', 'A2') ;
    
    xlswrite(XlsFilesave, RGBHeadings,'areas', 'A1') ;
    xlswrite(XlsFilesave,areas,'areas', 'A2') ;
    
    GeneralSetting = [numberOfBlobs, redThresholdLow , redThresholdHigh , greenThresholdLow , greenThresholdHigh,  blueThresholdLow ,blueThresholdHigh];
    Headings = {'numberOfBlobs','redThresholdLow' , 'redThresholdHigh' , 'greenThresholdLow' , 'greenThresholdHigh',  'blueThresholdLow' ,'blueThresholdHigh'};
    xlswrite(XlsFilesave, GeneralSetting,'GeneralSetting', 'A2') ;
    xlswrite(XlsFilesave, Headings,'GeneralSetting', 'A1') ;
    xlswrite(XlsFilesave, {'Procents of the image coverd by the chosen color:'},'GeneralSetting', 'A4') ;
    xlswrite(XlsFilesave, Area_Procent,'GeneralSetting', 'A5') ;
    
    xlswrite(XlsFilesave, {'Minimale area [pixel] of the blobs'},'GeneralSetting', 'A7') ;
    xlswrite(XlsFilesave, answer2,'GeneralSetting', 'A8') ;
end
%% Displaying the end of the computation
disp('*************************************************************')
disp('*****              Analysis has finnished               *****')
disp('*****          Check the xls file for results           *****')
disp('*************************************************************')
