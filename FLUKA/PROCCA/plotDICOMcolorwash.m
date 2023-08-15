function [varargout] = plotDICOMcolorwash(varargin)

%varargout

imgBG = varargin{1};
imgFG = varargin{2};
width = varargin{3};
start = varargin{4};

%% Remove the extra dimension
imgBG = squeeze(imgBG)';
imgFG = squeeze(imgFG)';

%% Pad the image and calculate the start offsets
s = max(size(imgBG') .* width);
offset = ((size(imgBG') .* width) - s) / 2;

%% Interpolate to center imgBackground, padding with zeros
imgBG = interp2(imgBG, (width(1):width(1):s)/width(1) + ...
    ((size(imgBG,2) * width(1)) - s)/(2 * width(1)), ...
    (width(2):width(2):s)'/width(2) + ...
    ((size(imgBG,1) * width(2)) - s)/(2 * width(2)), 'spline', 0);

%% Interpolate to center imageA, padding with zeros
imgFG = interp2(imgFG, (width(1):width(1):s)/width(1) + ...
    ((size(imgFG,2) * width(1)) - s)/(2 * width(1)), ...
    (width(2):width(2):s)'/width(2) + ...
    ((size(imgFG,1) * width(2)) - s)/(2 * width(2)), 'spline', 0);

%% Create spatial reference object based on the start and width inputs
RI = imref2d(size(imgBG), [start(1) + offset(1) start(1) + ...
    offset(1) + size(imgBG,2) * width(1)], [-(start(2) - offset(2) + ...
    size(imgBG,1) * width(2)) -start(2) + offset(2)]);

%% Plot 2D image data

% doseMC_bw = imgFG .* (imgFG > 0);
% alpha = (imgFG > 0.1) .* 0.6; % EGENTLIG
% 10^(floor(log10(min(min(imgFG)))))

alpha = (imgFG > 0.0001) .* 0.6; % 0.001 for proton and 0.000001 for photon
imgFG = imgaussfilt(imgFG, [3 3]); % [3 3]

% figure();
% imshow(imgFG, ...
%     'DisplayRange', [0 10^(floor(log10(max(max(imgFG)))))], ...
%     'ColorMap', colormap(jet(4096)));

h = 1;

if nargin > 4
    
    dose_ref            = varargin{5};
    imgFG               = 100 .* (imgFG ./ dose_ref); % max(max(imgFG)) ;
    %imgFG(imgFG > 100)  = 100;
    
    disp(['1: ' num2str(max(max(imgFG)))])
    disp(['2: ' num2str(10^(floor(log10(max(max(imgFG)))))) ])
    disp(['3: ' num2str(ceil(max(max(imgFG))/10)*10)])

    % h = 1;
    h = figure();
%     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(ind2rgb(gray2ind((imgBG) / 2048, 256), colormap('gray')), RI)
    hold on
    ih = imshow(imgFG, RI, ...
        'DisplayRange', [0 10^(floor(log10(max(max(imgFG)))))], ...
        'ColorMap', colormap(jet(4096)));
    set(ih, 'AlphaData', alpha);
    cb = colorbar;
    cb.Label.String = 'Dose (%)'; % 'Dose (nGy per primary)'; % 
    cb.FontSize = 30;
    caxis([0 100]) % 400 10^(floor(log10(max(max(imgFG)))))
    xlabel('cm', 'FontSize', 30)
    ylabel('cm', 'FontSize', 30)
    % set(gca, 'FontSize', 30)
    
%     h2 = figure();
%     set(h2, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%     imshow(ind2rgb(gray2ind((imgBG) / 2048, 256), colormap('gray')), RI)
%     
%     h1 = figure();
%     set(h1, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%     imshow(imgFG, ...
%         'DisplayRange', [0 10^(floor(log10(max(max(imgFG)))))], ...
%         'ColorMap', colormap(jet(4096)));
%     roi = drawfreehand(); % ; drawcircle();
%     bw_temp = roi.createMask();
%     figure(); imshow(imgFG .* bw_temp, [])
%     roi_values = reshape(imgFG .* bw_temp, 1, []);
%     roi_values(roi_values == 0) = [];
%     roi_values(isnan(roi_values)) = [];
%     [mean(roi_values) std(roi_values) prctile(roi_values, [5 25 50 75 95])]
    
    
end

%% Store the return variables
if nargout >= 1; varargout{1} = h; end
if nargout >= 2; varargout{2} = imgFG; end
if nargout >= 3; varargout{3} = RI; end
if nargout >= 4
    
    ROI             = [340 350 210 220]; % [340 350 260 270];
    dose_ref        = mean2( imgFG(ROI(3):ROI(4), ROI(1):ROI(2)) );
    varargout{4}    = dose_ref;
    
%     figure();
%     imshow(imgFG(ROI(3):ROI(4), ROI(1):ROI(2)), ...
%         'DisplayRange', [0 10^(floor(log10(max(max(imgFG)))))], ...
%         'ColorMap', colormap(jet(4096)));

    imgFG = 100 .* (imgFG ./ dose_ref);
%     sd              = std2( imgFG(ROI(3):ROI(4), ROI(1):ROI(2)) )
%       imgFG(imgFG > 150)  = 150;
%       ceil(max(max(imgFG))/10)*10
    10^(floor(log10(max(max(imgFG)))))

    % h = 1;
    h = figure();
    % set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(ind2rgb(gray2ind((imgBG) / 2048, 256), colormap('gray')), RI)
    hold on
    rectangle('Position', [ROI(1) * RI.PixelExtentInWorldX, ...
       ROI(3) * RI.PixelExtentInWorldY, ...
       (ROI(2)-ROI(1)) * RI.PixelExtentInWorldX, ...
       (ROI(4)-ROI(3)) * RI.PixelExtentInWorldY], ...
       'EdgeColor', 'k', 'LineWidth', 2)
    ih = imshow(imgFG, RI, 'DisplayRange', [0 10^(floor(log10(max(max(imgFG)))))], ...
        'ColorMap', colormap(jet(4096)));
    set(ih, 'AlphaData', alpha);
    cb = colorbar;
    cb.Label.String = 'Dose (%)'; % 'Dose (nGy per primary)'; 
    caxis([0 10^(floor(log10(max(max(imgFG)))))]) % 400
    xlabel('cm')
    ylabel('cm')
    set(gca, 'FontSize', 20)
    hold off
   
end

end