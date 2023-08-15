function [tform] = registerImage(img_fixed, img_moving)

[optimizer, metric] = imregconfig('multimodal');

tform = imregtform(img_moving, img_fixed, 'rigid', optimizer, metric);

end
