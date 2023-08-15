function [] = plot3Dsurf(x, y, z, X, Y, Z1, Z2, str_y, str_lgd, path)

% Compute density
% Put points into 3D bins; xyzBinNum is an nx3 matrix containing
% the bin ID for n values in xyz for the [x,y,z] axes.
xbins     = linspace(min(x(:)), max(x(:)), 100+1);
ybins     = linspace(min(y(:)), max(y(:)), 19+1);
zbins     = linspace(min(z(:)), max(z(:)), 7+1); % 7
xyzBinNum = [discretize(x(:), xbins), discretize(y(:), ybins), discretize(z(:), zbins)];

% bin3D is a mx3 matrix of m unique 3D bins that appear
% in xyzBinNum, sorted.  binNum is a nx1 vector of bin
% numbers identifying the bin for each xyz point. For example,
% b=xyz(j,:) belongs to bins3D(b,:).
[bins3D, ~, binNum] = unique(xyzBinNum, 'rows');
% density is a mx1 vector of integers showing the number of
% xyz points in each of the bins3D. To see the number of points
% in bins3D(k,:), density(k).
density = histcounts(binNum,[1:size(bins3D,1), inf])';
% Compute bin centers
xbinCnt = xbins(2:end)-diff(xbins)/2;
ybinCnt = ybins(2:end)-diff(ybins)/2;
zbinCnt = zbins(2:end)-diff(zbins)/2;

if Z2 ~= 0
    h = figure();
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    bubblechart3(xbinCnt(bins3D(:,1)), ybinCnt(bins3D(:,2)), ...
        zbinCnt(bins3D(:,3)), density, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'red', 'MarkerFaceAlpha', .4)
    % plot3(x, y, z, 'r.')
    % scatter3(x, y, z, 'ro', 'MarkerSize', 15)
    hold on
    surf(X, Y, Z1, 'FaceColor', 'green', 'EdgeColor', 'k') % HER
    surf(X, Y, Z2, 'FaceColor', 'blue', 'EdgeColor', 'k', 'FaceAlpha', 0.5)
    hold off
    xlabel('Dose (Gy)')
    ylabel(str_y)
    zlabel('Surviving Colony Count (per quadrat)')
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    zlim([min(z) max(z)])
    grid on; box on; camlight; view(40,35)
    legend({['Colony survival data $(N=' num2str(length(z)) ')$'], ...
        'LQ fit $(\exp(\lambda_0 -\alpha D - \beta D^2))$', str_lgd}, ...
        'Interpreter', 'latex', 'Location', 'NorthEast')
    bubblelegend('Quadrat count', 'Location', 'NorthEast')
    set(gca, 'FontSize', 16)
    % saveas(h, path)
else
    h = figure();
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    bubblechart3(xbinCnt(bins3D(:,1)), ybinCnt(bins3D(:,2)), ...
        zbinCnt(bins3D(:,3)), density, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'red', 'MarkerFaceAlpha', .4)
    % plot3(x, y, z, 'r.')
    % scatter3(x, y, z, 'ro', 'MarkerSize', 15)
    hold on
    surf(X, Y, Z1, 'FaceColor', 'green', 'EdgeColor', 'k') % HER
    hold off
    xlabel('Dose (Gy)')
    ylabel(str_y)
    zlabel('Surviving Colony Count (per quadrat)')
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    zlim([min(z) max(z)])
    grid on; box on; camlight; view(40,35)
    legend({['Colony survival data $(N=' num2str(length(z)) ')$'], ...
        'LQ fit $(\exp(\lambda_0 -\alpha D - \beta D^2))$'}, ...
        'Interpreter', 'latex', 'Location', 'NorthEast')
    bubblelegend('Quadrat count', 'Location', 'NorthEast')
    set(gca, 'FontSize', 16)
    % saveas(h, path)
end

end