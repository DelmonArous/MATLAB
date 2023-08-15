dx = linspace(-10,10,60);
dy = linspace(-10,10,60);

[x,y] = meshgrid(dx,dy);

u = cos(x).*sin(y);
v = -sin(x).*cos(y);
quiver(x,y,u,v,1.5)
axis equal