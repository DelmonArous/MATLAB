function [struct] = selectPrincipalComponentChannel(Z, rows, cols, dims, ...
    bw, tform)

% Specify pixel spatial relationships, i.e. distance between the pixel of 
% interest and its neighbor
offsets = [0 1; -1 1; -1 0; -1 -1];
% offsets = [0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1; 1 0; 1 1];
Rfixed = imref2d(size(bw));

for i = 1:dims
    
    struct.(sprintf('channel%i', i)).img = ...
        double(reshape(Z(:,i), rows, cols));
    
    img = double(reshape(Z(:,i), rows, cols));
    
    %     % Opening-by-Reconstruction
    %     img_e       = imerode(img, strel('disk', 25));
    %     img_obr     = imreconstruct(img_e, img);
    %
    %     % Opening-Closing-by-Reconstruction
    %     img_obrd    = imdilate(img_obr, strel('disk', 25));
    %     img_obrcbr  = imreconstruct(imcomplement(img_obrd), ...
    %         imcomplement(img_obr));
    %     bckg        = imcomplement(img_obrcbr);
    %
    %     % Suppress background
    %     img = img - bckg;
    %
    %     % Gaussian filtering
    %     img = imgaussfilt(img, [2 2]);
        
%     bw(bw == 0) = NaN;
%     img = img .* bw;
    img = (img - min(min(img))) ./ abs(max(max(img)) - min(min(img)));
    img = imwarp(img, tform, 'OutputView', Rfixed);
    img_CLAHE = adapthisteq(img, 'NumTiles', [16 16], 'ClipLimit', 0.008); % , 'Distribution', 'rayleigh'); 0.004
%     bw(bw == 0) = NaN;
    img_CLAHE = img_CLAHE .* bw;
    
%     img_diff = img - img_CLAHE;
%     img_diff = (img_diff - min(min(img_diff))) ./ abs(max(max(img_diff)) - min(min(img_diff)));
%     figure(); imshow(img_CLAHE, [])

    img_entropyfilt = entropyfilt(img);

    img = img_CLAHE;
    
    % Estimate metrics
    struct.(sprintf('channel%i', i)).rawEntropy    = entropy(img);
    struct.(sprintf('channel%i', i)).meanEntropy   = mean2(img_entropyfilt);
    struct.(sprintf('channel%i', i)).medianEntropy = median(img_entropyfilt, 'all');
    
    [glcm, SI] = graycomatrix(img, 'Offset', offsets, 'NumLevels', 128, ...
        'GrayLimits', []); %  minthresh maxthresh
    glcmstats                                    = graycoprops(glcm);
    struct.(sprintf('channel%i', i)).Entropy     = entropy(rescale(SI));
    struct.(sprintf('channel%i', i)).Contrast    = mean(glcmstats.Contrast);
    struct.(sprintf('channel%i', i)).Correlation = mean(glcmstats.Correlation);
    struct.(sprintf('channel%i', i)).Energy      = mean(glcmstats.Energy);
    struct.(sprintf('channel%i', i)).Homogeneity = mean(glcmstats.Homogeneity);
    
end

%% Channel selection by minimizing GLCM contrast value

% Minimum GLCM contrast
if struct.channel1.Contrast < struct.channel2.Contrast
    struct.channelopt.img       = struct.channel1.img;
    %         struct.channelopt.Contrast  = struct.channel1.GLCM.Contrast;
    struct.channelopt.channel   = 'PCA1';
else
    struct.channelopt.img       = struct.channel2.img;
    %         struct.channelopt.Contrast  = struct.channel2.GLCM.Contrast;
    struct.channelopt.channel   = 'PCA2';
end
    
% Maximum GLCM energy
% if strcmp(str, 'PCA')
%
%     if struct.channel1.Energy > struct.channel2.Energy
%         struct.channelopt.img       = struct.channel1.img;
% %         struct.channelopt.Contrast  = struct.channel1.GLCM.Contrast;
%         channel                     = '1';
%     else
%         struct.channelopt.img       = struct.channel2.img;
% %         struct.channelopt.Contrast  = struct.channel2.GLCM.Contrast;
%         channel                     = '2';
%     end
%
% else
%
%     if struct.channel1.Energy > struct.channel2.Energy && ...
%             struct.channel1.Energy > struct.channel3.Energy
%         struct.channelopt.img       = struct.channel1.img;
% %         struct.channelopt.Contrast  = struct.channel1.GLCM.Contrast;
%         channel                     = '1';
%     elseif struct.channel2.Energy > struct.channel1.Energy && ...
%             struct.channel2.Energy > struct.channel3.Energy
%         struct.channelopt.img       = struct.channel2.img;
% %         struct.channelopt.Contrast  = struct.channel2.GLCM.Contrast;
%         channel                     = '2';
%     elseif struct.channel3.Energy > struct.channel1.Energy && ...
%             struct.channel3.Energy > struct.channel2.Energy
%         struct.channelopt.img       = struct.channel3.img;
% %         struct.channelopt.Contrast  = struct.channel3.GLCM.Contrast;
%         channel                     = '3';
%     end
%
% end

% Minimum entropy
% if strcmp(str, 'PCA')
%
%     if struct.channel1.Entropy < struct.channel2.Entropy
%         struct.channelopt.img       = struct.channel1.img;
% %         struct.channelopt.Contrast  = struct.channel1.GLCM.Contrast;
%         channel                     = '1';
%     else
%         struct.channelopt.img       = struct.channel2.img;
% %         struct.channelopt.Contrast  = struct.channel2.GLCM.Contrast;
%         channel                     = '2';
%     end
%
% else
%
%     if struct.channel1.Entropy < struct.channel2.Entropy && ...
%             struct.channel1.Entropy < struct.channel3.Entropy
%         struct.channelopt.img       = struct.channel1.img;
% %         struct.channelopt.Contrast  = struct.channel1.GLCM.Contrast;
%         channel                     = '1';
%     elseif struct.channel2.Entropy < struct.channel1.Entropy && ...
%             struct.channel2.Entropy < struct.channel3.Entropy
%         struct.channelopt.img       = struct.channel2.img;
% %         struct.channelopt.Contrast  = struct.channel2.GLCM.Contrast;
%         channel                     = '2';
%     elseif struct.channel3.Entropy < struct.channel1.Entropy && ...
%             struct.channel3.Entropy < struct.channel2.Entropy
%         struct.channelopt.img       = struct.channel3.img;
% %         struct.channelopt.Contrast  = struct.channel3.GLCM.Contrast;
%         channel                     = '3';
%     end
%
% end

end