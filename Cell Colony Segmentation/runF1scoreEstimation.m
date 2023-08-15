clc
clear all
close all

%% Define directories for images and segmentation results
path_img        = 'C:\Users\delmo\Desktop\Jacob\Test_v4\Test Images';
path_bw_seg     = 'C:\Users\delmo\Desktop\Jacob\Test_v4\Test Colony Segmentation Masks';
% path_segdata    = 'C:\Users\delmo\Desktop\Jacob\Test\Test Colony Segmentation Data';
path_GT         = 'C:\Users\delmo\Desktop\Jacob\Test_v4\GT';
destpath        = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Manuscripts\EBT3 GRID Paper\Figures\New Figures';

%% Store file paths for each directory
filelist_img        = getAllFiles(path_img);
filelist_bw_seg     = getAllFiles(path_bw_seg);
folderlist_GT       = getAllFolders(path_GT);

%% Region of interest (ROI) size
ROI_size_mm = [8 8]; % in mm
dpi         = 1200;
px_size     = 25.4/dpi; % mm/px

%%

struct = {};
counter = 0;

for j = 1:length(folderlist_GT)

    counter = counter + 1;
    filelist_GT = getAllFiles(folderlist_GT{j});

    for i = 1:length(filelist_img)

        [~, fn_img, ~]      = fileparts(filelist_img{i});
        [~, fn_bw_seg, ~]   = fileparts(filelist_bw_seg{i});
        [~, fn_GT, ~]       = fileparts(filelist_GT{i});

%         strcmp(fn_img, fn_GT) && strcmp(fn_img, fn_bw_seg(1:end-8))

        if strcmp(fn_img, fn_GT) && strcmp(fn_img, fn_bw_seg(1:end-8))

            [struct.observer{counter}.img{i}.img, ...
                struct.observer{counter}.img{i}.bw_seg, ...
                struct.observer{counter}.img{i}.ROI_img, ...
                struct.observer{counter}.img{i}.ROI_bw_seg, ...
                struct.observer{counter}.img{i}.ROI_xy_GT, ...
                struct.observer{counter}.img{i}.ROI_bb, ...
                struct.observer{counter}.img{i}.N_seg, ...
                struct.observer{counter}.img{i}.N_GT, ...
                struct.observer{counter}.img{i}.precision, ...
                struct.observer{counter}.img{i}.recall, ...
                struct.observer{counter}.img{i}.F1, ...
                struct.observer{counter}.img{i}.filename] = ...
                computeF1score(filelist_img{i}, filelist_bw_seg{i}, ...
                filelist_GT{i}, dpi, ROI_size_mm);

        end

    end

end

%% Plot results

colorspec   = {'r*', 'g*', 'b*'};
dose        = [2 10 0 2 5 5 5 10 10 0]; % [2 2 5 10 0 0 2 10 0 2 5 5 5 10 10 0]; %  
p           = {};
precision   = [];
recall      = [];
F1          = [];
N_seg       = [];
N_GT        = [];

for i = 1:length(filelist_img)

    pvec        = [];
    lgdstr      = {'Segmentation outline', ...
        'Observer 1', 'Observer 2', 'Observer 3'};

    % Plot ROI superimposed on colony image
    h_img = figure(100*i);
    set(h_img, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(struct.observer{counter}.img{i}.img, [])
    hold on
    visboundaries(struct.observer{counter}.img{i}.bw_seg, ...
        'LineStyle', '-', 'Color', 'black', 'LineWidth', 1)
    rectangle('Position', struct.observer{counter}.img{i}.ROI_bb, ...
        'LineStyle', '--', 'EdgeColor', 'red', 'LineWidth', 4.5)
    errorbar(0.1*size(struct.observer{counter}.img{i}.img, 2), ...
        0.95*size(struct.observer{counter}.img{i}.img, 1), ...
        round(1/(px_size)), 'horizontal', 'k.', 'LineWidth', 1.5, ...
        'CapSize', 6);
    text(0.1*size(struct.observer{counter}.img{i}.img, 2), ...
        0.95*size(struct.observer{counter}.img{i}.img, 1), {'2 mm'}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'Fontsize', 12);
    hold off
%     set(gca, 'FontSize', 14)
%     saveas(h_img, fullfile(destpath, sprintf('%s_img.png', ...
%                         struct.observer{counter}.img{i}.filename)))

    % Plot ROI image 
    h_ROI = figure(i);
    set(h_ROI, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(struct.observer{counter}.img{i}.ROI_img)
    hold on

    boundaries = bwboundaries(struct.observer{counter}.img{i}.ROI_bw_seg);
    for k = 1:length(boundaries)
        B = boundaries{k};
        p{1} = plot(B(:,2), B(:,1), 'Color', 'black', 'LineWidth', 1);
    end
    pvec = [pvec p{1}];

    for j = 1:length(struct.observer)

        p{j+1} = plot(struct.observer{j}.img{i}.ROI_xy_GT(:,1), ...
            struct.observer{j}.img{i}.ROI_xy_GT(:,2), colorspec{j}, ...
            'MarkerSize', 18);
%         p{j+1}.MarkerFaceAlpha = 0.2;
        pvec = [pvec p{j+1}];

        N_seg(i,j)      = struct.observer{j}.img{i}.N_seg;
        N_GT(i,j)       = struct.observer{j}.img{i}.N_GT;
        precision(i,j)  = struct.observer{j}.img{i}.precision;
        recall(i,j)     = struct.observer{j}.img{i}.recall;
        F1(i,j)         = struct.observer{j}.img{i}.F1;

        lgdstr{j+1}     = lgdstr{j+1}; % [lgdstr{j+1} ...
%             ', precision=' num2str(round(precision(i,j),2)) ...
%             ', recall=' num2str(round(recall(i,j),2))]; 
        
    end

    hold off
    lgd = legend(pvec, lgdstr, 'Interpreter', 'latex');
    title(lgd, ['$D=' num2str(dose(i)) '$ Gy'], ...
        'Interpreter', 'latex', 'FontWeight', 'bold')
%     title(lgd, ['$D=' num2str(dose(i)) '$ Gy, $F_1=' ...
%         num2str(round(mean(F1(i,:)),2)) '\pm' ...
%         num2str(round(std(F1(i,:)),2)) '$'], ...
%         'Interpreter', 'latex', 'FontWeight', 'bold')
    set(gca, 'FontSize', 28)
%     saveas(h_ROI, fullfile(destpath, sprintf('%s_ROI.png', ...
%                         struct.observer{counter}.img{i}.filename)))

end

%%
path_img        = 'C:\Users\delmo\Desktop\Jacob\Temporary\Images';
path_bw_seg     = 'C:\Users\delmo\Desktop\Jacob\Temporary\Segmentation masks';
path_flask      = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB.tif';

filelist_img        = getAllFiles(path_img);
filelist_bw_seg     = getAllFiles(path_bw_seg);

f = uifigure;

for i = 1:length(filelist_img)

%     close all

    bw_seg = logical(readmatrix(filelist_bw_seg{i}));
    [img_in, ~] = imread(filelist_img{i});
    img_original = img_in(:,:,1:3);

    graychannel     = double(rgb2gray(img_in(:,:,1:3))); 
    img_gray        = ( graychannel - min(graychannel(:)) ) ./ ...
        abs( max(graychannel(:)) - min(graychannel(:)) );
    img_in          = double(img_gray);

    img_norm2 = normalizationMinMax(img_in);

    if mean(img_norm2(:)) > 0.45
        img_norm = normalizationMinMax(imcomplement(double(img_in(:,:,1))));
        img_norm(img_norm > 0.9) = 0.9;
    else
        img_norm = img_norm2;
    end

    
    [~, bw_border, tform] = getCellContainerBorder(path_flask, img_norm, f);
    for j = 1:3
        img_original(:,:,j) = imwarp(img_original(:,:,j), tform, ...
            'OutputView', imref2d(size(img_norm)));
%         img_original(:,:,j) = imrotate(img_original(:,:,j), 90);
    end
%     bw_seg = imrotate(bw_seg, 90);

    % Plot
    h_img = figure(1000*i);
    set(h_img, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(img_original, [])
    hold on
    visboundaries(bw_seg, 'LineStyle', '-', 'Color', 'black', 'LineWidth', 1)
    errorbar(0.1*size(img_original, 2), 0.95*size(img_original, 1), ...
        round(1/(px_size)), 'horizontal', 'k.', 'LineWidth', 1.5, ...
        'CapSize', 6);
    text(0.1*size(img_original, 2), 0.95*size(img_original, 1), ...
        {'2 mm'}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'Fontsize', 12);
    hold off

end

%%
deltax      = 0.5;
deltay      = 0.25;
ratio       = mean(N_GT,2)./mean(N_seg,2);
stdratio    = ratio .* sqrt( (std(N_GT,[],2) ./ mean(N_GT,2)).^2 + ...
    (std(N_seg,[],2) ./ mean(N_seg,2)).^2 );
precision_temp      = sum(precision, 2);
recall_temp         = sum(recall, 2);
F1_temp             = sum(F1, 2);
precision_std_temp  = std(precision, [], 2);
recall_std_temp     = std(recall, [], 2);
F1_std_temp         = std(F1, [], 2);

tbl = table(dose', ratio, 'VariableNames', {'D', 'f'});
mdl = fitlm(tbl, 'f ~ D + D^2', 'Intercept', true);
dosefit     = linspace(min(dose), max(dose), 1001);
ratiofit    = mdl.Coefficients.Estimate(1) + ... 
    mdl.Coefficients.Estimate(2) .* dosefit + ...
    mdl.Coefficients.Estimate(3) .* dosefit.^2;

% fittype2 = fittype('a + b*x + c*x^2');
% [fitobject2, gof2]  = fit(dose', ratio, fittype2);
% coeffval            = coeffvalues(fitobject2);
% coeffconfint        = confint(fitobject2);

figure();
h1 = errorbar(dose, ratio, stdratio, 'bo', 'MarkerSize', 6);
% plot(fitobject1, dose, ratio)
hold on
h2 = plot(dosefit, ratiofit, 'k-');
xlabel('Dose (Gy)', 'Interpreter', 'latex', 'FontWeight','bold')
ylabel('$\overline{N}_{GT} / N_{seg}$', 'Interpreter', 'latex', 'FontWeight','bold')
xlim([0-deltax max(dose)+deltax])
% ylim([0-deltay max(mean(N_GT,2)./mean(N_seg,2))+deltay])
legend([h1, h2], 'Data (\mu \pm \sigma)', 'Fit (\it{u + vD + wD^2})')
text(3.5, 1.7, '\bf{Coefficients (with standard deviations)}', ...
    'Interpreter', 'latex', 'FontWeight','bold')
text(4.6, 1.5, ['$u=' num2str(round(mdl.Coefficients.Estimate(1),2)) ...
    '\pm' num2str(round(mdl.Coefficients.SE(1),3)) '$'], ...
    'Interpreter', 'latex', 'FontWeight','bold')
text(4.6, 1.3, ['$v=' num2str(round(mdl.Coefficients.Estimate(2),2)) ...
    '\pm' num2str(round(mdl.Coefficients.SE(2),3)) '$'], ...
    'Interpreter', 'latex', 'FontWeight','bold')
text(4.6, 1.1, ['$w=' num2str(round(mdl.Coefficients.Estimate(3),3)) ...
    '\pm' num2str(round(mdl.Coefficients.SE(3),4)) '$'], ...
    'Interpreter', 'latex', 'FontWeight','bold')
text(4.6, 0.9, ['$R^2=' num2str(round(mdl.Rsquared.Ordinary,3)) ...
    ', \,\, RMSE=' num2str(round(mdl.RMSE,3)) '$'], ...
    'Interpreter', 'latex', 'FontWeight','bold')
grid on
set(gca, 'FontSize', 16)
hold off


% figure();
% h1 = errorbar(dose, F1_avg, F1_std, 'ro', 'MarkerSize', 6);
% hold on
% h2 = errorbar(dose, precision_avg, precision_std, 'go', 'MarkerSize', 6);
% h3 = errorbar(dose, recall_avg, recall_std, 'bo', 'MarkerSize', 6);
% xlabel('Dose (Gy)', 'Interpreter', 'latex', 'FontWeight','bold')
% ylabel('Performance metric', 'Interpreter', 'latex', 'FontWeight','bold')
% xlim([0-deltax max(dose)+deltax])
% legend([h1, h2, h3], '$F_1$ score', 'Precision', 'Recall', ...
%     'Location', 'SouthWest', 'Interpreter', 'latex', 'FontWeight', 'bold')
% set(gca, 'FontSize', 16)
% hold off


i = 0;
for doselevel = unique(dose)

    i = i + 1;
    ind_temp        = find(dose == doselevel);

    precision_avg(i)    = sum(precision_temp(ind_temp)) / (3*length(ind_temp));
    recall_avg(i)       = sum(recall_temp(ind_temp)) / (3*length(ind_temp));
    F1_avg(i)           = sum(F1_temp(ind_temp)) / (3*length(ind_temp));
    precision_std(i)    = sqrt(sum(precision_std_temp(ind_temp).^2) / length(ind_temp));
    recall_std(i)       = sqrt(sum(recall_std_temp(ind_temp).^2) / length(ind_temp));
    F1_std(i)           = sqrt(sum(F1_std_temp(ind_temp).^2) / length(ind_temp));

end


% text(10, 800, ...
%   sprintf('\\it{a = %0.2g (%0.2g, %0.2g) (CI)\n  b = %0.2g (%0.2g, %0.2g) (CI)\n }', ...
%   c, cint(2)-c));

% figure()
% yyaxis left
% % errorbar(dose-0.07, mean(F1,2), std(F1,[],2), 'o', 'MarkerSize', 6)
% plot(dose, mean(F1,2), 'o', 'MarkerSize', 6)
% xlabel('Dose (Gy)', 'Interpreter', 'latex', 'FontWeight','bold')
% ylabel('$F_1$ score', 'Interpreter', 'latex', 'FontWeight','bold')
% xlim([0-deltax max(dose)+deltax])
% ylim([0-deltay 1+deltay])
% 
% yyaxis right
% hold on
% % errorbar(dose+0.07, mean(N_GT,2)./mean(N_seg,2), stdratio, 'o', 'MarkerSize', 6)
% plot(dose, mean(N_GT,2)./mean(N_seg,2), 's', 'MarkerSize', 6)
% hold off
% ylabel('$N_{GT} / N_{seg}$', 'Interpreter', 'latex', 'FontWeight','bold')
% ylim([0-deltay max(mean(N_GT,2)./mean(N_seg,2))+deltay])
% set(gca, 'FontSize', 12)

% %% Plot results
% figure();
% imshow(ROI_img{1}{1})
% % title(fn)
% hold on
% visboundaries(ROI_seg{1}{1}, 'Color', 'red', 'LineWidth', 1)
% plot(xy_GT{1}{1}(:,1), xy_GT{1}{1}(:,2), 'b*')
% hold off
