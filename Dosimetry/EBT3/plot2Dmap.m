function [] = plot2Dmap(imgs, px_size, range, label)

% RA = imref2d(size(imgs,1,2), px_size, px_size);
temp_img = sum(imgs, 3) ./ size(imgs, 3);
temp_img = imrotate(temp_img, 90);
RA = imref2d(size(temp_img,1,2), px_size, px_size);

h = figure();
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
imshow(temp_img, RA, range, 'colormap', jet(4096)) % 5 Gy dose delivered
shading interp
c = colorbar;
c.Label.String = label;
xlabel('cm'); ylabel('cm');
% title(lgdstr)
grid on
set(gca, 'FontSize', 17)

end