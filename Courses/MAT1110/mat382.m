dx = -10:0.5:10;
dy = -10:0.5:10;

[x,y] = meshgrid(dx, dy);

u = -y./(x.^2 + y.^2);
v = x./(x.^2 + y.^2);

hold on
quiver(x, y, u, v, 1.5)
streamline(x,y, u,v, 1,0)
hold off
