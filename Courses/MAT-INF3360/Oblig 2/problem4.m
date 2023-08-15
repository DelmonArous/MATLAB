close all;
clear all;
clc;

Tmax = 1.;
nx = 99;                    % number of space steps  
nt = 9999;                  % number of time steps
deltax = 1./(nx+1);         % deltax = x(2) - x(1)
deltat = Tmax/(nt+1);

x = linspace(0, 1, nx+2);
t = linspace(0, Tmax, nt+2);
v_expl = zeros(nx+2, nt+2);

f = @(x)(x.*(1-x));
g = @(t)(0.);
a = @(x)(1.);
u = @(x,t)(f(x-a(x).*t).*(x>a(x).*t) + ...
        g(t-x./a(x)).*(x<=a(x).*t));   % Analytical solution

% Initial condition
for j = 1:nx+2
    x_j = (j-1)*deltax;
    v_expl(j,1) = f(x_j);
end

% Boundary condition
for m = 2:nt+2
    t_m = (m-1)*deltat;
    v_expl(1,m) = g(t_m);
end

% Implementation of the explicit method
for m = 1:nt+1
    for j = 2:nx+2
        r_j = a((j-1)*deltax)*(deltat/deltax);
        v_expl(j,m+1) = r_j*v_expl(j-1,m) + (1 - r_j)*v_expl(j,m); 
    end
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
S2 = surf(X, T, v_expl');
set(S2, 'EdgeColor', 'none', 'FaceColor', 'interp');
view([-60 25]);
set(gcf, 'Color', 'w', 'Units', 'pixels', 'Position', [200 200 700 600]);
xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);
zlabel('u(x,t)', 'FontSize', 12);

figure();
contourf(X, T, abs(U-v_expl'));
colorbar();
xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);