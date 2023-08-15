function [I] = iterativeFeedbackOtsu(img, otsuvec, a_min, a_max)

% Initialize
I   = false(size(img));

%% Combination of masks generated for each Otsu value 
for i = 1:length(otsuvec)
    
    bw = im2bw(img, otsuvec(i));
    bw = bwpropfilt(bw, 'Area', [a_min/2 2*a_max]);
    bw = imfill(bw, 'holes');
    bw = bwpropfilt(bw, 'Area', [a_min/2 2*a_max]);
    bw = imfill(bw, 'holes');
    bw = imdilate(bw, ones(3));
    
%     if parameters.areavec(2) < 200
%         bw = imopen(bw, strel('disk', ceil(1/parameters.otsuvector(i))));
%     else
%         bw = imopen(bw, strel('disk', ceil(3/parameters.otsuvector(i))));
%     end
    
    bw = bwpropfilt(bw, 'Eccentricity', [0 0.98]);
    bw = bwpropfilt(bw, 'Area', [a_min/2 2*a_max]);
    I = I | bw;
    
end

end