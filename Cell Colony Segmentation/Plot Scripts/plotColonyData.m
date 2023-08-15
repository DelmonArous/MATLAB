function [] = plotColonyData(area, red, green, blue, gray, pca1, pca2, ...
    destpath, filename)

area = area.';
red = red.'; green = green.'; blue = blue.'; gray = gray.';
pca1 = pca1.'; pca2 = pca2.';

% %% Plot RGB colony values  
% figure();
% scatter3(red, green, blue, 25, area, 'filled', 'MarkerEdgeColor', 'black')
% xlabel('Red'), ylabel('Green'), zlabel('Blue');
% xlim([0 1]), ylim([0 1]), zlim([0 1]);
% axis equal
% colormap(jet(4096))
% cb = colorbar;
% cb.Label.String = 'Colony Area (mm^2)';
% caxis([0 ceil(max(area(:)))])
% set(gca, 'FontSize', 14)

%% Inter-colony correlation and plot PCA1 and PCA2 colony values
[r, pval] = corr(pca2, pca1, 'type', 'Pearson');
str = strcat('n = ', num2str(length(area)), ...
    ', R^2 = ', num2str(round(r^2,3)), ...
    ', p-value = ', num2str(round(pval,7)));
X = [ones(length(pca1), 1) pca1];
beta = X\pca2;

h = figure();
hold on
scatter(pca1, pca2, 25, area, 'filled', 'MarkerEdgeColor', 'black')
plot(pca1, X*beta, 'r-');
hold off
xlabel('1st Principal Component'), ylabel('2nd Principal Component'); 
colormap(jet(4096))
cb = colorbar;
cb.Label.String = 'Colony Area (mm^2)';
caxis([0 ceil(max(area(:)))])
shading interp
lgd = legend('Mean colony value', 'Location', 'best');
title(lgd, str)
set(gca, 'FontSize', 14)
saveas(h, fullfile(destpath, sprintf('%s-CorrPCA.fig', filename)))

%% Plot variable correlation
X = [area red green blue gray pca1 pca2];
varnames = {'Area', 'Red', 'Green', 'Blue', 'Gray', 'PCA1', 'PCA2'};
h = figure();
[R, Pvalue] = corrplot(X, 'testR', 'on', 'varNames', varnames);
% title(['Pearson''s Correlation Matrix, n = ' ...
%     num2str(length(area)) ' colonies'])
saveas(h, fullfile(destpath, sprintf('%s-Corr.fig', filename)))


end