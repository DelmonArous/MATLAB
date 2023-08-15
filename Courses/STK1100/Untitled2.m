clear all

t = 8.0;
ant_alpha = 0:11;
ant_periode = [57 203 383 525 532 408 273 139 45 27 10 6];
periode_tid = t.*ant_periode;
P = [];

lambda1 = t*sum((ant_alpha./ant_periode))/length(ant_alpha)

lambda = 3.8715;

rel_frek = ant_periode./sum(ant_periode);

for k = 1:12
    P(k) = exp(-lambda)*(lambda)^(k-1)/factorial(k-1);
end

plot(ant_alpha, P, 'b')
hold on
plot(ant_alpha, rel_frek, 'r')
ylabel('Sannsynlighet')
xlabel('Antall \alpha-partikler')
legend('Poissonfordeling', 'Relativ frekvens')
