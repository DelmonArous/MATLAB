clear all; clf; 
% Physical constants
F = 400;     % N 
m = 80;      % kg
rho = 1.293; % kg/m^3
A = 0.45;    % m^2
C_D = 1.2;
w = 0;       % m/s
dt = 0.01;   % s
time = 9;    % s
t_c = 0.67;  % s
f_c = 488;   % N
f_v = 25.8;  % sN/m
% Initial conditions
v0 = 0;
x0 = 0;
% Numerical initialization
n = time/dt;
x = zeros(n, 1);
v = zeros(n, 1);
a = zeros(n, 1);
t = zeros(n, 1);
% Set initial values
x(1) = x0;
v(1) = v0;
% Integration loop
for i = 1:n-1
    F_C = f_c*exp(-(t(i)/t_c)^2);
    F_V = f_v*v(i);
    A = A*(1-0.25*exp(-(t(i)/t_c)^2));
    D = 0.5*A*rho*C_D*(v(i)-w)^2;
    F_net = F + F_C - F_V - D;
    a(i+1) = F_net/m;
    v(i+1) = v(i) + a(i+1)*dt;
    x(i+1) = x(i) + v(i+1)*dt;
    t(i+1) = t(i) + dt; 
end
% Plot results
subplot(3,1,1);
plot(t,x);
text(8.68, 100.1052,'\leftarrow 100m','HorizontalAlignment','left')
xlabel('t [s]');
ylabel('x [m]');
subplot(3,1,2);
plot(t,v);
xlabel('t [s]');
ylabel('v [m/s]');
subplot(3,1,3);
plot(t,a);
xlabel('t [s]');
ylabel('a [m/s^2]');