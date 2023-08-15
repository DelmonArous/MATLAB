clear all; clf; 
% Physical constants
F = 400;     % N 
m = 80;      % kg
rho = 1.293; % kg/m^3
A = 0.45;    % m^2
C_D = 1.2;
w = 0;       % m/s
dt = 0.01;   % s
time = 50;   % s
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
    a(i+1) = (F - 0.5*rho*C_D*A*(v(i))^2)/m;
    v(i+1) = v(i) + a(i+1)*dt;
    x(i+1) = x(i) + v(i+1)*dt;
    t(i+1) = t(i) + dt;
end
disp(max(v)) % 33.8492 m/s is the numerical terminal velocity 
% Plot results
subplot(3,1,1);
plot(t,x);
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