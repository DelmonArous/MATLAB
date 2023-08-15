m = 0.1;    % kg
k = 100;    % N/m
b = 0.1;    % m
g = 9.8;    % m/s^2
time = 2;   % s
dt = 0.01;  % s
u = 0.0;    % m/s
omega = sqrt(k/m);

n = ceil(time/dt);

v0 = 0.1;   % m/s
x0 = 0.0;   % m

t = zeros(n,1);
x = zeros(n,1);
v = zeros(n,1);
x_exact = zeros(n,1);

v(1) = v0;
x(1) = x0;

for i = 1:n-1
    x_b = u*t(i) + x0 + b;
    a = (k/m)*(x_b-x(i)-b);
    v(i+1) = v(i) + a*dt;
    x(i+1) = x(i) + v(i+1)*dt;
    x_exact(i+1) = (v0/omega)*sin(omega*t(i));
    t(i+1) = t(i) + dt;
end

plot(t, x, '-', t, x_exact, '--')
xlabel('t [s]')
ylabel('x [m]')
title('Plot of the motion on the block, both numerical and exact motion')
legend('Numerical','Exact')