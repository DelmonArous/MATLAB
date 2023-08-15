[x,y,u,v] = velfield(51);

quiver(x,y,u,v,1.5)
xlabel('x-akse')
ylabel('y-akse')
title('Hastighetsfelt: v = (cos(x)sin(y), -sin(x)cos(y)')
axis equal