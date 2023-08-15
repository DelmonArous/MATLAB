close all;
clear all;
clc;

dybde = [0.5 1 1.5 2 3 5 10 15];
foton_6MV = (5/3)*[1.092 1.12 1.21 1.2 1.15 1.05 0.8 0.6];
foton_15MV = [1.25 1.31 1.63 1.73 1.78 1.68 1.35 1.08];

plot(dybde, foton_6MV, '-rx', dybde, foton_15MV, '-bo');
grid on
xlabel('Vanndybde [cm]');
ylabel('Indusert strøm [nA]');
legend('6 MV fotoner', '15 MV fotoner');
axis([0 16 0 2.2])

figure();
horisontal_avstand = [0 4 4.5 5 5.5 6 7];
indusert_strom = [1 0.963 0.897 0.645 0.29 0.14 0.075];
plot(horisontal_avstand, indusert_strom, '-ro');
grid on
xlabel('Horisontal forflytning fra sentrum [cm]');
ylabel('Relativ indusert strøm');
axis([0 8 0 1.2])