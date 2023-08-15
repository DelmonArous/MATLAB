function [img_norm] = MinMaxNorm(img, param)

% Background check
if median(median(img)) < max(max(img))
    npx_white = numel(find(img >= median(median(img))));
    npx_black = numel(find(img < median(median(img))));
else
    npx_white = numel(find(img >= mean(mean(img))));
    npx_black = numel(find(img < mean(mean(img))));
end

if strcmp(param.pca, 'PCA2')
    px_flag = npx_white < npx_black;
else
    px_flag = npx_white > npx_black;
end

% if ~isempty(npx_white) && ~isempty(npx_black) && (npx_white > npx_black)
%     flag = true;
% else
%     flag = false;
% end
 
% img_vec                 = double(img(:));
% maxthresh               = prctile(img_vec, 99);
% minthresh               = prctile(img_vec, 1);
% img(img > maxthresh)    = maxthresh;
% img(img < minthresh)    = minthresh;

if ~isempty(npx_white) && ~isempty(npx_black) && px_flag % flag == true
    img         = imcomplement(img);
    img_norm    = im2double(img);
    img_norm    = (img_norm - min(min(img_norm))) / ...
        abs(max(max(img_norm)) - min(min(img_norm)));
else
    img_norm = (img - min(min(img)))/abs(max(max(img)) - min(min(img)));
end

end