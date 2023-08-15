dx = -10:0.1:10;
dy = -10:0.1:10;
[x,y] = meshgrid(dx,dy);

for i = 1:length(x)
    if x(i) == 0 && y(i) == 0
        z = 0;
    else
        z = (x.^3.*y-x.*y.^3)./(x.^2 + y.^2);
    end
end

meshc(x,y,z)