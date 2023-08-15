clear all;
clc;

Tmax = 1.;
l1 = -2.;
l2 = 3.;
nx = 25;                    % number of space steps  
nt = 5000;                  % number of time steps
deltax = (l2-l1)/(nx+1);    % deltax = x(2) - x(1)
deltat = Tmax/(nt+1);
r = deltat/deltax^2

x = linspace(l1, l2, nx+2);
t = linspace(0, Tmax, nt+2);

v_expl = zeros(nx+2, nt+2);
v_impl = zeros(nx+2, nt+2);
u_array = zeros(nx+2, nt+2);
q_array = zeros(nx+2, nt+2);

u = @(x,t)(exp(t)+x);       % Analytical solution
q = @(x,t)(11*exp(t)+10*x);

% Initial condition
f = @(x)(1. + x);
for j = 2:nx+1
    x_j = (j-1)*deltax;
    v_expl(j,1) = f(x_j);
    v_impl(j,1) = f(x_j);
end

% Boundary condition
a = @(t)(exp(t) - 2.);
b = @(t)(exp(t) + 3.);
for m = 1:nt+2
    t_m = (m-1)*deltat;
    v_expl(1,m) = a(t_m);
    v_expl(nx+2,m) = b(t_m);
    v_impl(1,m) = a(t_m);
    v_impl(nx+2,m) = b(t_m); 
end

% Implementation of the explicit method
for m = 1:nt+1
    for j = 2:nx+1
        v_expl(j,m+1) = 4*r*v_expl(j-1,m) + ...
            (1-8*r-10*deltat)*v_expl(j,m) + 4*r*v_expl(j+1,m) ...
            + deltat*q((j-1)*deltax,(m-1)*deltat);
    end
end

% Implementation of the implicit method
alpha(1:nx-1) = -4*r; beta(1:nx) = 1+8*r+10*deltat; gamma(1:nx-1) = -4*r; 
A = diag(alpha,-1) + diag(beta,0) + diag(gamma,1);
for m = 2:nt+2
    t_m = (m-1)*deltat;
    v_prev = v_impl(2:nx+1, m-1) + deltat*q(x(2:nx+1),t_m).';
    v_impl(2:nx+1, m) = A\v_prev;
end

plot(x, v_expl(:,500), '--', x, v_impl(:, 500), '.-', x, u(x,t(500)), ':')
title(['Analytical and numerical solutions at t = ' t(500)])
xlabel('x')
ylabel('u(x,t)')
%legend('Explicit','Implicit','Analytical')