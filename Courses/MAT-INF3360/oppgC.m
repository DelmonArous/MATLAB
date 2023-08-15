clear all;
clc;

n = [10, 20, 40, 80, 160];
h = 1./(n+1.);

I1_approx = zeros(length(n),1);
I2_approx = zeros(length(n),1);
I1_exact = 1/6.;
I2_exact = sqrt(2)/3.;
alpha1 = zeros(length(n)-1,1);
alpha2 = zeros(length(n)-1,1);

F = @(x)(x^5);
G = @(x)(sqrt(abs(x-0.5)));
    
for i = 1:length(n)
    I1_approx(i) = trapz(F, 0, 1, n(i));
    I2_approx(i) = trapz(G, 0, 1, n(i));
end

format long
I1_approx, I2_approx

abserror1 = abs(I1_exact - I1_approx)
abserror2 = abs(I2_exact - I2_approx)
relerror1 = abserror1./I1_exact;
relerror2 = abserror2./I2_exact;

for i = 1:(length(n)-1)
    alpha1(i) = log( abserror1(i)/abserror1(i+1) )/log( h(i)/h(i+1) );
    alpha2(i) = log( abserror2(i)/abserror2(i+1) )/log( h(i)/h(i+1) );
end

alpha1, alpha2