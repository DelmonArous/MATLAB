function [datamap, RI] = read2DMap(filename)

%% Read 2D dose map
fileID = fopen(filename,'r');
formatSpec = '%f %f %f %f';
dataarray_size = [4 Inf];   
dataarray = fscanf(fileID, formatSpec, dataarray_size);
dataarray =  dataarray';

%% Store each column into vectors
x = dataarray(:,1); y = dataarray(:,2); % in cm
data = dataarray(:,3);

%% Rearrange the data vector into a 2D map
xi = unique(x); yi = unique(y);
[X,Y] = meshgrid(xi, yi);
datamap = reshape(data, size(X));
% datamap = datamap .* 1.602176462*10^(-7) * 10^9; % in nGy
%dosemap = flip(dosemap, 2);

% Get rid of (dosemap <= 0), (dosemap == Inf) and (dosemap == NaN)
datamap(datamap <= 0) = 0; 
datamap(datamap == Inf) = 0;
datamap(isnan(datamap)) = 0;

% Find min/max 
% dose_min = min(datamap(:));
% dose_max = max(datamap(:));
% dose_range = dose_max - dose_min;

% Convert dose data to log space
% dosemap = log10(dosemap);

% Normalize the log scaled dosemap so that the colorbar scale will be correct
% dosemap = dose_min + (dose_range) .* (dosemap - log10(dose_min)) ./ ...
%     (log10(dose_max/dose_min)); % self scale data
% L = [0.00000001 0.0000001 0.000001 0.00001 0.0001 0.001 0.01 0.1 1];
% l = dose_min + (dose_range) .* (log10(L) - log10(dose_min)) ./ ...
%     (log10(dose_max/dose_min)); % Tick mark positions

% Create a spatial referencing object associated with the image, and
% use the referencing object to set the x- and y-axes limits in the
% world coordinate system
% conversion = ( abs(min(xi(:))) + abs(max(xi(:))) ) / (length(xi)-1); % in mm/pixel
% sizex = size(dosemap, 2);
% sizey = size(dosemap, 1);
% xmax = sizex * conversion;
% ymax = sizey * conversion;
RI = imref2d(size(datamap));
RI.XWorldLimits = [min(xi(:)) max(xi(:))];
RI.YWorldLimits = [min(yi(:)) max(yi(:))];

%% Plot
% h = figure();
% set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% imshow(dosemap, RI, [], 'colormap', jet(256)) % [min(dosemap(:)) max(dosemap(:))] 
% c = colorbar;
% c.Label.String = 'Dose (nGy/primary)';
% % set(c, 'Ytick', l, 'YTicklabel', L);
% set(gca, 'FontSize', 16)
% xlabel('cm'); ylabel('cm');
% shading interp

end
