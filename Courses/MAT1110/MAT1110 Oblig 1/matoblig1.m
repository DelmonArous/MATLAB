r = -2:0.1:2;
s = -1.4:0.1:1.4;
t = -1.4:0.1:1.4;

n = zeros(length(t), 1);

r_x = t;
r_y = t.^2;
r_z = t.^3;

[x, y] = meshgrid(r, s);
z = y.^2./x;

hold on
plot3(r_x, r_y, r_z)
plot3(n, n, z)
surf(x, y, z)
hold off

axis([-2 2 -2 2 -4 4]);
xlabel('x-akse');
ylabel('y-akse');
zlabel('z-akse');