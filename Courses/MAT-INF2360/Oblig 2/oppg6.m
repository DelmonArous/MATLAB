close all;
clear all;
clc;

omega = linspace(0,2*pi,1000);
lambda_g0 = (1./8)*(cos(2*omega) + 4*cos(omega) + 3.);
lambda_g1 = (1./64)*(5*cos(5*omega) + 20*cos(4*omega) ...
    + cos(3*omega) - 96*cos(2*omega) - 70*cos(omega) + 140);

figure();
plot(omega, abs(lambda_g0), omega, abs(lambda_g1))
xlabel('\omega')
ylabel('|\lambda(\omega)|')
legend(['\lambda_{G_0}'], ['\lambda_{G_1}'])