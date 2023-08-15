function [I_norm] = img_normalization(img)

img = double(img);
prc_vec = [1 99];
I_low = prc_vec(1); I_high = prc_vec(2);
type = 'minmax';
pixrev  = false;

if isempty(pixrev)
    if median(median(im)) < max(max(im))
        npx_white = numel(find(im >= median(median(im))));
        npx_black = numel(find(im < median(median(im))));
    else
        npx_white = numel(find(im >= mean(mean(im))));
        npx_black = numel(find(im < mean(mean(im))));
    end
    
    if  (npx_white > npx_black) && ~isempty(npx_white) && ~isempty(npx_black)
        pixrev = true;
    end  
else
    pixrev = true;
end

if pixrev == true
    img = imcomplement(img);
    I_norm = im2double(img);
    I_norm = (img - min(min(img)))/abs(max(max(img)) - min(min(img)));
else
    I_norm = (img - min(min(img)))/abs(max(max(img)) - min(min(img)));
end


end