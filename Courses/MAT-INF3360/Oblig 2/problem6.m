close all;
clear all;
clc;

Tmax = 1.;
nx = 999;                    % number of space steps  
nt = 99999;                  % number of time steps
deltax = 1./(nx+1)        % deltax = x(2) - x(1)
deltat = Tmax/(nt+1)

x = linspace(0, 1, nx+2);
t = linspace(0, Tmax, nt+2);
v_impl = zeros(nx+2, nt+2);

f = @(x)(x.*(1-x));
g = @(t)(0.);
a = @(x)(1.);
u = @(x,t)(f(x-a(x).*t).*(x>a(x).*t) + ...
        g(t-x./a(x)).*(x<=a(x).*t));   % Analytical solution

% Initial condition
for j = 1:nx+2
    x_j = (j-1)*deltax;
    v_impl(j,1) = f(x_j);
end

% Boundary condition
for m = 2:nt+2
    t_m = (m-1)*deltat;
    v_impl(1,m) = g(t_m);
end

% Implementation of the implicit method
alpha(1:nx) = -deltat/deltax; beta(1:nx+1) = 1 + deltat/deltax;
A = diag(alpha,-1) + diag(beta,0);
for m = 1:nt+1
    t_m = (m-1)*deltat;
    v_prev = v_impl(2:nx+2, m);
    v_impl(2:nx+2, m+1) = A\v_prev;
end

[X,T] = meshgrid(x,t);
U = u(X,T);

figure();
S1 = surf(X, T, U);
set(S1, 'EdgeColor', 'none', 'FaceColor', 'interp');
view([-60 25]);
set(gcf, 'Color', 'w', 'Units', 'pixels', 'Position', [200 200 700 600]);
xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);
zlabel('u(x,t)', 'FontSize', 12)

figure();
S2 = surf(X, T, v_impl');
set(S2, 'EdgeColor', 'none', 'FaceColor', 'interp');
view([-60 25]);
set(gcf, 'Color', 'w', 'Units', 'pixels', 'Position', [200 200 700 600]);
xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);
zlabel('u(x,t)', 'FontSize', 12);

figure();
contourf(X, T, abs(U-v_impl'));
colorbar();
xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);