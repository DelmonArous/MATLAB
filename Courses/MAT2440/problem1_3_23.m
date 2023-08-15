[X Y] = meshgrid(0:0.1:3, -1:0.1:4);
dY = X.^2 + Y.^2 -  1;
dX = ones(size(dY));
L = sqrt(dX.^2+dY.^2);
quiver(X,Y, dX./L, dY./L), axis tight