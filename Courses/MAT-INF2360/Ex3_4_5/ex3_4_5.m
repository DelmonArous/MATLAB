clear all;
close all;
clc;

c = 0.1;
omega = linspace(0, 2*pi, 1000);
lambda_k4 = 1 + c.*exp(-1i*omega) +(c^2).*exp(-2i*omega)...
    + (c^3).*exp(-3i*omega) + (c^4).*exp(-1i*omega);

lambda_k8 = 1 + c.*exp(-1i*omega) +(c^2).*exp(-2i*omega)...
    + (c^3).*exp(-3i*omega) + (c^4).*exp(-1i*omega)...
    + (c^5).*exp(-5i*omega) + (c^6).*exp(-6i*omega)...
    + (c^7).*exp(-7i*omega) + (c^8).*exp(-8i*omega);

figure;
plot(omega, lambda_k4, omega, lambda_k8);
legend('k=4', 'k=8')
