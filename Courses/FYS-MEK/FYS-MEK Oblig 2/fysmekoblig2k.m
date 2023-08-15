clear all; clf; 
% Physical constants
m = 80;      % kg
rho = 1.293; % kg/m^3
A = 0.45;    % m^2
C_D = 1.2;
w = 0;       % m/s
dt = 0.01;   % s
time = 10;   % s
t_c = 0.67;  % s
f_c = 488;   % N
f_v = 25.8;  % sN/m
% Initial conditions
v0 = 0;
% Numerical initialization
n = time/dt;
v = zeros(n, 1);
t = zeros(n, 1);
F = zeros(n, 1);
F_C = zeros(n, 1);
F_V = zeros(n, 1);
D = zeros(n, 1);
% Set initial values
v(1) = v0;
% Integration loop
for i = 1:n-1
    F(i) = 400;
    F_C(i) = f_c*exp(-(t(i)/t_c)^2);
    F_V(i) = f_v*v(i);
    D(i) = 0.5*A*(1-0.25*exp(-(t(i)/t_c)^2))*rho*C_D*(v(i)-w)^2;
    F_net = F(i) + F_C(i) - F_V(i) - D(i);
    a(i+1) = F_net/m;
    v(i+1) = v(i) + a(i+1)*dt;
    t(i+1) = t(i) + dt; 
end
% Plot results
plot(t, F,'-', t, F_C,'--', t, F_V,':', t, D,'-.');
legend('F','F_C','F_V','D')
xlabel('t [s]');
ylabel('F [N]');

% The driving force F is constant throughout the motion, it is the only
% force that is independent.
% From the net-force-model we see that the initial driving forc is 
% decreasing exponentialy with the time t, so its effect on the runner 
% will therfore fall.
% However, the driving force F_V from the model will grow over time 
% since the velocity of the runner enhances over time. This driving 
% force will then stabilize when the runner reaches the terminal velocity,
% hence the effect of this driving force will be at it's peek when the
% terminal velocity is reached.
% The model for the air resistance D is velocity-dependent, thus it's
% effect on the runner will become approximately constant when the change
% of the velocity becomes close to zero(i.e. when the runner reaches
% the terminal velocity).