clear all;
close all;
clc;

delta_min = [4.748*10^4 58983 1.446*10^5];
delta_max = [7.506*10^4 95392 2.5616*10^5];
dose = [10 20 68];

plot(dose, delta_min, '*', dose, delta_max, 'x');
xlabel('Dose [Gy]');
ylabel('Intensitet');
p1 = polyfit(dose, delta_min, 1);
yfit1 = polyval(p1, dose);
p2 = polyfit(dose, delta_max, 1);
yfit2 = polyval(p2, dose);
hold on

plot(dose, yfit1, '-r', dose, yfit2, '-g');

