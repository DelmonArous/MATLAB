% Modify from here -->
N = 3;          % nr of balls
m = 0.05;       % kg
k = 10000;      % N/m
q = 4;
d = 0.015;      % m
v0 = 0.5;       % m/s
time = 0.7;     % s
dt = 0.001;     % s
% Base variables
n = ceil(time/dt);
x = zeros(n,N);
v = zeros(n,N);
F = zeros(n,N);
t = zeros(n,1);
% Initial conditions, equally spaced
for j = 1:N
    x(1,j) = d*(j-1);
end
v(1,1) = v0;
for i = 1:n-1
  % Find force on each block, j
  % First, force from block to the left
  for j = 2:N
    dx = x(i,j) - x(i,j-1);
    F(i,j) = F(i,j) + force(dx,d,k,q);
  end
  % Second, force from block to the right
  for j = 1:N-1
    dx = x(i,j+1) - x(i,j);
    F(i,j) = F(i,j) - force(dx,d,k,q);
  end
  % Euler-Cromer step
  for j = 1:N
    a = F(i,j)/m;
    v(i+1,j) = v(i,j) + a*dt;
    x(i+1,j) = x(i,j) + v(i+1,j)*dt;
  end
  %
  % The Euler-Cromer step above can also be vectorized
  % through the following implementation (which is faster)
  % Euler-Cromer vectorized step
  % a = F(i,:)/m;
  % v(i+1,:) = v(i,:) + a*dt;
  % x(i+1,:) = x(i,:) + v(i,:)*dt;
  %
  t(i+1) = t(i) + dt;
end
% Plot results
for j = 1:N
  plot(t,v(:,j));
  if j==1
    hold on;
  end
  if j==N
    hold off
  end
end
xlabel('t [s]');
ylabel('v [m/s]');
title(['Plot of the velocities on each ball with physical constants:'... 
    'N=',num2str(N),' balls, m=',num2str(m),'kg, k=',num2str(k),...
    'N/m, q=',num2str(q),', d=',num2str(d),'m, v_0=',num2str(v0),...
    'm/s, dt=',num2str(dt),'s, time=',num2str(time),'s'])

