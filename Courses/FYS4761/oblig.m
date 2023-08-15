close all;
clear all;
clc;

mu_60keV = [0.0 0.395 0.260 0.227 0.196 0.181 0.164];
CT_tall = [-950.63 913.26 335.40 120.61 -32.93 -87.47 -173.22];

p1 = polyfit(mu_60keV, CT_tall, 1);
y1 = polyval(p1, mu_60keV);

figure;
plot(mu_60keV, CT_tall, 'o', mu_60keV, y1, '-');
xlabel('\mu [cm^{-2}]');
ylabel('CT-tall [HU]');
legend('Data','Tilpasset kurve')