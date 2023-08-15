close all;
clear all;
clc;

Tmax = 1.;
nx = 49;                    % number of space steps  
nt = 4999;                  % number of time steps
deltax = 1./(nx+1)        % deltax = x(2) - x(1)
deltat = Tmax/(nt+1)
s = deltat/deltax^2;
epsilon = 0.1;

x = linspace(0, 1, nx+2);
t = linspace(0, Tmax, nt+2);
v = zeros(nx+2, nt+2);

f = @(x)(x.*(1-x));
g = @(t)(0.);
h = @(t)(0.);
a = @(x)(1.);

% Initial condition
for j = 1:nx+2
    x_j = (j-1)*deltax;
    v(j,1) = f(x_j);
end

% Boundary condition
for m = 2:nt+2
    t_m = (m-1)*deltat;
    v(1,m) = g(t_m);
    v(nx+2,m) = h(t_m);
end

% Implementation of the scheme
alpha(1:nx-1) = -0.5*deltat/deltax - epsilon*s; 
beta(1:nx) = 1 + 2*epsilon*s;
gamma(1:nx-1) = 0.5*deltat/deltax - epsilon*s;
A = diag(alpha,-1) + diag(beta,0) + diag(gamma,1);
for m = 1:nt+1
    v_prev = v(2:nx+1, m);
    v(2:nx+1, m+1) = A\v_prev;
end

% Surface plot
[X,T] = meshgrid(x,t);

figure();
S = surf(X, T, v');
set(S, 'EdgeColor', 'none', 'FaceColor', 'interp');
view([-20 30]);
set(gcf, 'Color', 'w', 'Units', 'pixels', 'Position', [200 200 700 600]);
xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);
zlabel('u^\epsilon(x,t)', 'FontSize', 12);

