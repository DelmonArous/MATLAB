function [PV] = statsModel(img_model, bw, n_x, n_y, dxdy)

% Partition survival image into quadrats
temp_img_model_quadrat = mat2quadrat(img_model, dxdy);

% Initialization
img_model_quadrat = zeros(n_x, n_y);

for i = 1:n_x
    for j = 1:n_y
        img_model_quadrat(i,j) = mean(temp_img_model_quadrat{i,j}(:));
    end
end

% % Compute deviation maps
% img_diff_LQ_quadrat  = img_model_quadrat  - img_count_quadrat;
% img_diff_MLQ_quadrat = img_SF_MLQ_quadrat - img_count_quadrat;

% Get all relevant quadrat colony count data
PV = img_model_quadrat(bw); % regionprops(bw, img_model_quadrat, 'PixelValues');

end