k = 200;        % N/m
g = 9.81;       % m/s^2
tetha = pi/6;   % angle
m = 0.1;        % kg
L0 = 1.0;       % m

r0 = [L0*sin(tetha) -L0*cos(tetha)];
v0 = [0.0 0.0];

time = 10.0;
dt = 0.001;
n = time/dt;

r = zeros(n,2);
v = zeros(n,2);

r(1,:) = r0;
v(1,:) = v0;

for i = 1:n-1
    rr = norm(r(i,:));
    a = [0.0 -g] - (k/m)*(rr - L0)*(r(i,:)/rr);        
    v(i+1,:) = v(i,:) + dt*a;
    r(i+1,:) = r(i,:) + dt*v(i+1,:);
end

plot(r(:,1), r(:,2))
xlabel('x [m]')
ylabel('y [m]')
title(['Plot of the balls´ motion for the first 10s with dt=',num2str(dt)])