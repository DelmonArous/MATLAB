% Constants
g = 9.81;   % m/s^2
m = 0.1;    % kg
h = 0.001;   % 
time = 1.0; % s
dt = 0.001; % timestep
% Initial conditions
v0 = 0;
z0 = 1.0;
% Numerical initialization
n = time/dt;
z = zeros(n,1);
v = zeros(n,1);
a = zeros(n,1);
t = zeros(n,1);
% Set initial values
z(1) = z0;
v(1) = v0;
% Integration loop
for i = 1:n-1
    F = -m*g;
    
    a(i+1) = F/m;
    v(i+1) = v(i) + dt*a(i+1);
    z(i+1) = z(i) + dt*v(i+1);
    if (z(i)>h)&&(z(i+1)<h)
        v(i+1) = -v(i+1);
    end
    t(i+1) = t(i) + dt;
end

plot(v,z)
title('Plot of a bouncing ball')
xlabel('v [m/s]')
ylabel('z [m]')