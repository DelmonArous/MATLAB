clear all;
clc;

n = [8, 16, 32, 64, 128, 256, 1024];
h = 1./(n+1);

E_h1 = zeros(1, length(n));
E_h2 = zeros(1, length(n));
E_h3 = zeros(1, length(n));
E_h4 = zeros(1, length(n));
E_h5 = zeros(1, length(n));

alpha1 = zeros(1, length(n)-1);
alpha2 = zeros(1, length(n)-1);
alpha3 = zeros(1, length(n)-1);
alpha4 = zeros(1, length(n)-1);
alpha5 = zeros(1, length(n)-1);

f1 = @(x)(1.);
f2 = @(x)(x);
f3 = @(x)(x^2);
f4 = @(x)(exp(x));
f5 = @(x)(cos(x));

for i = 1:length(n)
    x = linspace(0, 1, n(i)+2);
    u1_exact = (x./2.).*(1-x);
    u2_exact = (x./6.).*(1-x.^2);
    u3_exact = (x./12.).*(1-x.^3);
    u4_exact = (exp(1)-1).*x - exp(x) + 1;
    u5_exact = (1-cos(1)).*x + cos(x) - 1;
    
    u1_approx = PoissonDirichletSolver(f1, n(i));
    u2_approx = PoissonDirichletSolver(f2, n(i));
    u3_approx = PoissonDirichletSolver(f3, n(i));
    u4_approx = PoissonDirichletSolver(f4, n(i));
    u5_approx = PoissonDirichletSolver(f5, n(i));
        
    E_h1(i) = max(abs(u1_exact-u1_approx));
    E_h2(i) = max(abs(u2_exact-u2_approx));
    E_h3(i) = max(abs(u3_exact-u3_approx));
    E_h4(i) = max(abs(u4_exact-u4_approx));
    E_h5(i) = max(abs(u5_exact-u5_approx));
end

format long
E_h1, E_h2, E_h3, E_h4, E_h5

for i = 1:(length(n)-1)
    alpha1(i) = log( E_h1(i)/E_h1(i+1) )/log( h(i)/h(i+1) );
    alpha2(i) = log( E_h2(i)/E_h2(i+1) )/log( h(i)/h(i+1) );
    alpha3(i) = log( E_h3(i)/E_h3(i+1) )/log( h(i)/h(i+1) );
    alpha4(i) = log( E_h4(i)/E_h4(i+1) )/log( h(i)/h(i+1) );
    alpha5(i) = log( E_h5(i)/E_h5(i+1) )/log( h(i)/h(i+1) );
end

alpha1, alpha2, alpha3, alpha4, alpha5

