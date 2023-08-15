function L = watershedBLOB(img, bw_BLOB, param, ws_vec, flag)

% The following variables are required for execution: 
%   img = input image that is either a grayscale or a single color (red) 
%       channel image containing a single inquired cell colony 
%       conglomeration (BLOB). 
%   bw_BLOB = binary mask containing the associated logical BLOB for 
%       watershed segmentation.
%   param = structure containing neccessary watershed parameters, such as
%       user-defined colony area and segmentation threshold
%   ws_vec = watershed vector containing the multi-threshold segmentation
%       steps.
%   flag = logical value to force a watershed segmentation procedure if 
%       input BLOB is still too big, regardless of shape (eccentricity and
%       circularity). 
%
% The following variables are returned upon successful completion when 
% input arguments are provided:
%   L = binary mask containing the segmented cell colonies.

warning('off', 'all')

%% Parameters
area = param.area;
medianarea = param.medianarea;
BLOBsize = nnz(bw_BLOB); % (find(bw_BLOB == 1));
boundary_BLOB = imdilate(bwperim(bw_BLOB), ones(3));

% Expected number of colonies to be found for each BLOB
E_m = ceil(BLOBsize/(medianarea));

% BLOB properties
stats_BLOB = regionprops(bw_BLOB, 'Area', 'Eccentricity', 'Circularity');
img_adapthisteq = adapthisteq(img);
img_adapthisteq_complement = imcomplement(img_adapthisteq);

% Initialization
L = cell(size(img,1), size(img,2), numel(ws_vec));
y = (-1) .* ones(length(ws_vec), 1);

% figure(); imshow(img_adapthisteq, [])
% figure(); imshow(bw_BLOB, [])

%% Watershed only if BLOBs are of more elongated shape

if flag || (stats_BLOB.Eccentricity > 0.4 || stats_BLOB.Circularity < 0.6)

    % Loop through each watershed segmentation threshold
    for i = 1:length(ws_vec)
        
%         img_adapthisteq(img_adapthisteq > ...
%             prctile(img_adapthisteq(:),95)) = ...
%             prctile(img_adapthisteq(:),95); % new
%         img_adapthisteq = imfill(img_adapthisteq, 'holes');
        
        % Detect only regional maxima that are larger than neighboring 
        % pixels by a selected intensity threshold
        bw_em = imextendedmax(img_adapthisteq, ws_vec(i));
        bw_em = imfill(bw_em, 'holes');
        bw_em = bwmorph(bw_em, 'thicken');

        % Remove all parts of the conglomerate BLOB that are located
        % outside the boundary of the BLOB
        bmask = boundary_BLOB & bw_em;
        bw_em = bw_em - bmask;
        bw_em = imclearborder(bw_em);

        % Euclidean distance transform
        D = -bwdist(~bw_em);
        bw_em = imextendedmin(D,2);
        D = imimposemin(D, bw_em);

        % Watershed
        Ld = watershed(D);
        bw_em = bw_BLOB;
        bw_em(Ld == 0) = 0;
        bw_em = ~bwperim(bw_em, 8);

        L{i} = logical(imclearborder(bw_em));
        L{i} = bwmorph(L{i}, 'clean');
%         L{i} = imopen(L{i}, strel('disk',1));
        if max(max(bwlabel(L{i}))) < 2
            L{i} = bw_BLOB;
        end

%         figure(); imshow(bwlabel(L{i}), [])

        % Count number of segmented colonies and store their features
        count = max(max(bwlabel(L{i})));
        stats = regionprops(logical(L{i}), 'Area', 'Eccentricity', ...
            'Circularity');
        
        % Initialize the quality vector
        mu_vec = (-1) .* ones(length(stats),1);
        
        % Calculate the quality criteria; mu = mu_1 * mu_2 * mu_3
        for j = 1:length(stats)
            mu_1   = evalmf(stats(j).Area, [area(1) area(2) ...
                max(area(3),area(2)*2) area(4)], 'pimf');
            mu_2 = evalmf(stats(j).Circularity, [0.35 0.6 0.9 1], 'pimf');
            mu_vec(j) = mu_1*mu_2;
        end
        
        if sum(sum(bw_BLOB)) > 2*area(4)
            mu_vec = median(mu_vec);
        else
            mu_vec = prod(mu_vec);
        end
        mu_3 = evalmf(count, [1 max(1,E_m) 2*E_m 3*E_m-1], 'pimf');
        y(i) = (mu_vec*mu_3);
       
    end
    
    % Select the result with maximum quality
    if max(y) > 0
        
        if E_m == 1
            ind  = find(y > 0);
        else
            ind  = find(y == max(y));
        end
        L  = L{ind(1)};
        
    else
        
        % Otherwise, a single-shot watershed segmentation is carried out 
        % by supressing outlier pixels on the higher side of the image 
        % intensity histogram
        img(img > prctile(img(:),95)) = prctile(img(:),75); % median(img(:));
        img = imfill(img, 'holes');
        bw_em = imextendedmax(img, 0.25*mean(ws_vec)); % 0.05
        bw_em = imfill(bw_em, 'holes');
        bw_em = bwmorph(bw_em, 'thicken');
        
        % Remove all parts of the conglomerate BLOB that are located
        % outside the boundary of the BLOB
        bmask = boundary_BLOB & bw_em;
        bw_em = bw_em - bmask;
        bw_em = imclearborder(bw_em);
        
        % Euclidean distance transform
        D = -bwdist(~bw_em);
        bw_em = imextendedmin(D,2);
        D = imimposemin(D, bw_em);
        
        % Watershed
        Ld = watershed(D);
        bw_em = bw_BLOB;
        bw_em(Ld == 0) = 0;
        bw_em = ~bwperim(bw_em, 8);
        L = logical(imclearborder(bw_em));
        L = bwmorph(L, 'clean');
        if max(max(bwlabel(L))) < 2
            L = bw_BLOB;
        end
        
    end
    
else
    
    L = zeros(size(img));
    
end

end