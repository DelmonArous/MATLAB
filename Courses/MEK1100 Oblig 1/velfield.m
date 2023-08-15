function [x,y,u,v] = velfield(n)

x = linspace(-10,10,n);
y = linspace(-10,10,n);

[x,y] = meshgrid(x,y);
u = cos(x).*sin(y);
v = -sin(x).*cos(y);

end