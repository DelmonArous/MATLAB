clear all;
clc;

Tmax = 1.;
nx = 43;                % number of space steps  
nt = 4000;              % number of time steps
deltax = 1./(nx+1);     % deltax = x(2) - x(1)
deltat = Tmax/(nt+1);
r = deltat/deltax^2

x = linspace(0, 1, nx+2);
t = linspace(0, Tmax, nt+2);
v_expl = zeros(nx+2, nt+2);
v_impl = zeros(nx+2, nt+2);

syms k;
u1 = @(x,t)(3.*exp(-pi^2*t)*sin(pi*x) + 5.*exp(-16.*pi^2*t)*sin(4*pi*x));
u2 = @(x,t)((4./pi)*symsum(exp(-((2*k-1)*pi)^2*t)*sin((2*k-1)*pi*x)/(2*k-1), 1,nx));
u3 = @(x,t)((2./pi)*symsum((-1)^(k+1)*exp(-(k*pi)^2*t)*sin(k*pi*x)/k, 1,nx));

% Initial condition
f1 = @(x)(3.*sin(pi*x) + 5.*sin(4*pi*x));   % example 3.1
f2 = @(x)(1.);                              % example 3.2
f3 = @(x)(x);                               % example 3.4
for j = 2:nx+1
    x_j = (j-1)*deltax;
    v_expl(j,1) = f2(x_j);
    v_impl(j,1) = f2(x_j);
end

% Boundary condition
for m = 1:nt+2
    v_expl(1,m) = 0.;
    v_expl(nx+2,m) = 0.;
    v_impl(1,m) = 0.;
    v_impl(nx+2,m) = 0.;    
end

% Implementation of the explicit method
for m = 1:nt+1
    for j = 2:nx+1
        v_expl(j,m+1) = r*v_expl(j-1,m) + (1.-2.*r)*v_expl(j,m) + r*v_expl(j+1,m);
    end
end

e = @(t)(max(abs(u1(x,t) - v_expl(:,t))));

% Implementation of the implicit method
a(1:nx-1) = -1.; b(1:nx) = 2.; c(1:nx-1) = -1.; 
A = (1./deltax^2)*(diag(a,-1) + diag(b,0) + diag(c,1));
for m = 1:nt+1
    v_prev = v_impl(2:nx+1, m);
    v_impl(2:nx+1, m+1) = (eye(nx)+deltat*A)\v_prev;
end

plot(x, v_expl(:,301), '--', x, v_impl(:,301), '.-', x, u2(x,0.1), ':')
title('Analytical and numerical solutions at t = 0.1')
xlabel('x')
ylabel('u(x,t)')
legend('Explicit','Implicit','Analytical')