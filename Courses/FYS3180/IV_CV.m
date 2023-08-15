clear all;
close all;
clc;

data_iv = load('iv_A_4E10.txt');
data_cv = load('cv_A_4E10.txt');

V_iv = data_iv(:,1);
I = data_iv(:,2).*10^(-12);

V_cv = data_cv(:,1);
C = data_cv(:,2)./10^(-12);

figure();
plot(V_iv, I, '-*')
xlabel('Bias voltage [V]','interpreter','latex','FontSize',13)
ylabel('Current [pA]','interpreter','latex','FontSize',13)

figure();
semilogy(V_iv, I, '-*')
xlabel('Bias voltage [V]','interpreter','latex','FontSize',13)
ylabel('$\log_{10}(I)$','interpreter','latex','FontSize',13)

figure();
plot(V_cv, C, '-*')
xlabel('Bias voltage [V]','interpreter','latex','FontSize',13)
ylabel('Capacitance [pF]','interpreter','latex','FontSize',13)