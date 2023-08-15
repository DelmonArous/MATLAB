dx = -10:0.1:10;
dy = -10:0.1:10;
[x,y] = meshgrid(dx,dy);
t = -10:0.1:10;

for i = 1:length(x)
    if x(i) == 0 && y(i) == 0
        z = 0;
    else
        z = x.^2.*y./x.^4 + y.^2;
    end
end

hold on
meshc(x,y,z)
plot3(t, t.^2, 0.5)
hold off