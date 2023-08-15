x = -5:0.5:5;
y = -5:0.5:5;

[x, y] = meshgrid(x,y);
vx = x.*y;
vy = y;

quiver(x, y, vx, vy, 1.5)
xlabel('x-akse')
ylabel('y-akse')

for i = -5:5
    streamline(x, y, vx, vy, i, -0.001)
end
for i = -5:5
    streamline(x, y, vx, vy, i, 0.001)
end