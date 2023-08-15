% clc
% clear all
% close all
% delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
% delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
% warning('off','all')


%% Directories

sourcepath_imgs         = 'C:\Users\delmo\Desktop\TestOpen\Images';
sourcepath_colonydata   = 'C:\Users\delmo\Desktop\TestOpen\ColonyData';
sourcepath_dish  = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB.tif';

%%

filelist_imgs   = getAllFiles(sourcepath_imgs);
dose            = repmat([0 2 5 10]', 1, 12)'; % dose in Gy
dose            = dose(:)';
f               = uifigure;

median_vec      = [];
mean_vec        = [];

for i = 1:length(filelist_imgs)
    
    [p, fn, ~] = fileparts(filelist_imgs{i});
    
    % If imread is successful, store the image information
    [img_in, ~]     = imread(filelist_imgs{i});
    img_gray        = double(rgb2gray(img_in)); 

    % Normalization
    [img_norm, ~] = normalizationMinMax(img_gray);
    
    bw_border = getCellContainerBorder(sourcepath_dish, img_norm, f);
     
    % Multiply the dish/flask mask with the colony image to extract out
    % solely free viable area in the container
    img_norm = bw_border .* img_norm;
    
%     figure(); imshow(img_norm, [])
%     figure(); imshow(bw_border, [])

    img_values = reshape(img_norm, 1, []);
    
    % Remove voxel values inside the mask of the structure,
    % in which are zero or NaN
    img_values(img_values == 0) = [];
    img_values(isnan(img_values)) = [];
    
    % Calculate the 1-99th R2* percentiles
    img_values_median   = median(img_values);
    img_values_mean     = mean(img_values);
    img_values_prctiles = prctile(img_values, 1:99);
    
    median_vec  = [median_vec img_values_median];
    mean_vec    = [mean_vec img_values_mean];
    
end

%%

filelist_colonydata = getAllFiles(sourcepath_colonydata);

N_colony_tmp = [2254 3242 3405 2171 3474 3505 2261 3484 3318 2175 3314 3512 ...
    1619 2932 2257 1658 3000 2863 1558 2931 2864 1691 2932 2837 647 ...
    1951 1729 658 2169 1721 659 1984 1864 682 2332 1877 2326 2987 2753 ...
    2052 3325 2541 1859 3204 3808 1726 3251 3031];
N_colony = [];

for i = 1:length(filelist_colonydata)

    data        = readtable(filelist_colonydata{i});
    N_colony    = [N_colony length(data.ColonyID)];

end

figure()
plot(dose, mean_vec, 'o')
xlabel('Dose (Gy)')
ylabel('mean(SI)') % 'Transmittance \it{T}' % '\it{OD}'
xlim([min(dose)-0.5 max(dose)+0.5])
grid on
set(gca, 'FontSize', 16)

%% Quadratic regression

eq_type = fittype('a*x^2 + b*x + c');
fit1    = fit(N_colony_tmp(1:end-12)', mean_vec(1:end-12)', eq_type);
[p, S]  = polyfit(N_colony_tmp(1:end-12)', mean_vec(1:end-12)', 3); 

x_fit = linspace(0, 4000, 100000000);
y_fit = polyval(p, x_fit);

N_colony_10Gy = round(interp1(y_fit, x_fit, mean_vec(end-11:end)));

% Plot
figure()
hold on
plot(N_colony_tmp(1:end-12), mean_vec(1:end-12), 'o')
plot(x_fit, y_fit, 'r-')
% plot(fit, '-')
hold off
xlabel('Colony Count')
ylabel('mean(SI)') % 'Transmittance \it{T}' % '\it{OD}'
grid on
set(gca, 'FontSize', 16)

figure()
hold on
plot(N_colony_tmp, mean_vec, 'o')
plot(x_fit, y_fit, 'r-')
% plot(fit, '-')
hold off
xlabel('Colony Count')
ylabel('mean(SI)') % 'Transmittance \it{T}' % '\it{OD}'
grid on
set(gca, 'FontSize', 16)

save('SIfile.mat', 'median_vec', 'mean_vec', 'N_colony_10Gy', 'filelist_imgs')