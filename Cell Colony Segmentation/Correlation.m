clear all
close all
clc

%% Dataset 1 (grayscale)

ACC_count       = [48 66 56 65 60 62 55 51 60 72 53 63 43 73 73 56]; 

ACC_pre = [0.71 0.71 0.73 0.89 0.85 0.73 0.62 0.59 0.67 0.78 0.62 0.63 0.77 0.67 0.66 0.77];
ACC_rec = [0.92 0.98 0.91 0.95 0.94 0.92 0.94 0.91 0.89 0.88 0.97 1.00 0.82 0.94 0.92 0.90];
ACC_F1  = [0.80 0.82 0.81 0.92 0.89 0.81 0.75 0.71 0.76 0.82 0.76 0.78 0.80 0.78 0.77 0.83];

[mean(ACC_pre) std(ACC_pre) min(ACC_pre) max(ACC_pre)]
[mean(ACC_rec) std(ACC_rec) min(ACC_rec) max(ACC_rec)]
[mean(ACC_F1) std(ACC_F1) min(ACC_F1) max(ACC_F1)]


%% Dataset 1

ACC_count       = [34 42 40 58 49 40 32 30 38 50 33 37 29 51 50 39];
ACS_count       = [30 45 45 61 56 46 34 33 34 42 36 35 29 46 47 41];
MCC_count       = [36.7 48.0 45.3 62.0 51.0 47.3 35.7 30.3 40.7 59.3 33.3 40.7 37.0 52.3 54.3 46.0];
MCC_count_SD 	= [2.4 1.3 1.8 2.7 2.0 2.2 2.9 1.8 1.8 3.6 1.1 2.2 2.7 1.8 1.8 2.0];
GT_count        = [37 48 45 61 54 49 36 33 45 64 34 40 40 52 52 48];

APE_ACC         = (abs(GT_count - ACC_count) ./ GT_count) .* 100;
APE_ACS         = (abs(GT_count - ACS_count) ./ GT_count) .* 100;
APE_MCC         = (abs(GT_count - MCC_count) ./ GT_count) .* 100;

[mean(APE_ACC) std(APE_ACC) min(APE_ACC) max(APE_ACC)]
[mean(APE_ACS) std(APE_ACS) min(APE_ACS) max(APE_ACS)]
[mean(APE_MCC) std(APE_MCC) min(APE_MCC) max(APE_MCC)]

% [r_ACC, pval_ACC] = corr(GT_count', ACC_count', 'type', 'Pearson');
% [r_ACS, pval_ACS] = corr(GT_count', ACS_count', 'type', 'Pearson');
% [r_MCC, pval_MCC] = corr(GT_count', MCC_count', 'type', 'Pearson');
% X = [ones(length(GT_count'), 1) GT_count'];
% beta_ACC = X\ACC_count';
% beta_ACS = X\ACS_count';
% beta_MCC = X\MCC_count';

% markerspec  = {'ro', 'go', 'bo', 'ko'};
% linespec    = {'r-', 'g-', 'b-', 'k-'};
% lgdstr      = {'Red channel', 'Green channel', 'Blue channel', ...
%     'Grayscale image'};

x = linspace(0, 100, 1000);

figure();
hold on
h1 = plot(GT_count, ACC_count, 'ro', 'MarkerSize', 12);
h2 = plot(GT_count, ACS_count, 'gs', 'MarkerSize', 12);
h3 = plot(GT_count, MCC_count, 'b^', 'MarkerSize', 12);
h4 = errorbar(GT_count, MCC_count, MCC_count_SD, 'b^', 'MarkerSize', 12);
h5 = plot (x, x, 'k--');
% h5 = plot(GT_count, X*beta_ACC, 'r-');
% h6 = plot(GT_count, X*beta_ACS, 'g-');
% h7 = plot(GT_count, X*beta_MCC, 'b-');
hold off
xlabel('GT count')
ylabel('Colony count') % 'Transmittance \it{T}'  '\it{netOD}' % '\it{OD}'
xlim([25 70])
ylim([25 70])
% lgd = legend([h1, h2, h4], ...
%     ['r = ' num2str(round(r_ACC,3)) ' (ACC)'], ...
%     ['r = ' num2str(round(r_ACS,3)) ' (AutoCellSeg)'], ...
%     ['r = ' num2str(round(r_MCC,3)) ' (MCC)'], ...
%     'Location', 'SouthEast');
lgd = legend([h1, h2, h4], 'ACC', 'AutoCellSeg', 'MCC', ...
    'Location', 'NorthWest');
title(lgd, 'Dataset 1')
grid on
set(gca, 'FontSize', 20)

%% Dataset 2

clear all
close all
clc

str         = {'\it{E. coli}', '\it{Klebs. pn.}', ...
    '\it{Pseud. ae.}', '\it{Staph. au.}'};
strainspec  = {'r', 'g', 'b', 'c'};
counterspec = {'o', 's'};
ACC_count 	= {[112 77 34], [67 46 28], [29 20 24], [12 103 85]};
ACS_count 	= {[114 86 40], [64 47 28], [29 20 24], [11 98 84]};
GT_count 	= {[116 80 32], [67 49 27], [29 20 25], [13 106 88]};
x           = linspace(0, 120, 1000);

figure();
hold on
for i = 1:length(str)
    
%     [r_ACC, pval_ACC] = corr(GT_count{i}', ACC_count{i}', 'type', 'Pearson');
%     [r_ACS, pval_ACS] = corr(GT_count{i}', ACS_count{i}', 'type', 'Pearson');

    h(i)                = plot(GT_count{i}, ACC_count{i}, ...
        [strainspec{i} counterspec{1}], 'MarkerSize', 12);
    h(length(str)+i)    = plot(GT_count{i}, ACS_count{i}, ...
        [strainspec{i} counterspec{2}], 'MarkerSize',12);
    lgdstr{i}               = ['ACC - ' str{i}];
    lgdstr{length(str)+i}   = ['AutoCellSeg - ' str{i}];
    
    APE_ACC         = (abs(GT_count{i} - ACC_count{i}) ./ GT_count{i}) .* 100;
    APE_ACS         = (abs(GT_count{i} - ACS_count{i}) ./ GT_count{i}) .* 100;
    
    [mean(APE_ACC) std(APE_ACC) min(APE_ACC) max(APE_ACC)]
    [mean(APE_ACS) std(APE_ACS) min(APE_ACS) max(APE_ACS)]
    
%     lgd = legend( ...
%         ['r = ' num2str(round(r_ACC,3)) ' (ACC)'], ...
%         ['r = ' num2str(round(r_ACS,3)) ' (AutoCellSeg)'], ...
%         'Location', 'SouthEast');
    
end
h(9) = plot (x, x, 'k--');
lgd = legend([h(1) h(2) h(3) h(4) h(5) h(6) h(7) h(8)], lgdstr, ...
    'Location', 'NorthWest');
xlabel('GT count')
ylabel('Colony count')
% xlim([min(GT_count{i})-5 max(GT_count{i})+5])
% ylim([min(GT_count{i})-5 max(GT_count{i})+5])
hold off
title(lgd, 'Dataset 2')
grid on
set(gca, 'FontSize', 20)

%% Dataset 2

% ACC_count 	= [112 77 34 67 46 28 29 20 24 12 103 85];
% ACS_count 	= [114 86 40 64 47 28 29 20 24 11 98 84];
% GT_count 	= [116 80 32 67 49 27 29 20 25 13 106 88];
% 
% x = linspace(0, 100, 1000);
% 
% % markerspec  = {'ro', 'go', 'bo', 'ko'};
% % linespec    = {'r-', 'g-', 'b-', 'k-'};

 
% figure();
% hold on
% h1 = plot(GT_count, ACC_count, 'ro');
% h2 = plot(GT_count, ACS_count, 'rs');
% h3 = plot(GT_count, MCC_count, 'r^');
% h4 = errorbar(GT_count, MCC_count, MCC_count_SD, 'r^');
% h5 = plot (x, x, 'k--');
% % h5 = plot(GT_count, X*beta_ACC, 'r-');
% % h6 = plot(GT_count, X*beta_ACS, 'g-');
% % h7 = plot(GT_count, X*beta_MCC, 'b-');
% hold off
% xlabel('GT count')
% ylabel('Automated count') % 'Transmittance \it{T}'  '\it{netOD}' % '\it{OD}'
% xlim([min(GT_count)-10 max(GT_count)+10])
% ylim([min(GT_count)-10 max(GT_count)+10])
% % lgd = legend([h1, h2, h4], ...
% %     ['r = ' num2str(round(r_ACC,3)) ' (ACC)'], ...
% %     ['r = ' num2str(round(r_ACS,3)) ' (AutoCellSeg)'], ...
% %     ['r = ' num2str(round(r_MCC,3)) ' (MCC)'], ...
% %     'Location', 'SouthEast');
% lgd = legend([h1, h2, h4], 'ACC', 'AutoCellSeg', 'MCC', ...
%     'Location', 'NorthWest');
% title(lgd, 'Dataset 1')
% grid on
% set(gca, 'FontSize', 16)
