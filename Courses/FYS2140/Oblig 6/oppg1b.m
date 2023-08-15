clear all;

x = linspace(-10,10,1000);

psi_0 = (1.0/pi)^(1/4)*exp((-1.0/2)*x.^2);
psi_1 = (1.0/pi)^(1/4)*sqrt(2.0)*x.*exp((-1.0/2.0)*x.^2);
psi_2 = (1.0/pi)^(1/4)*(1/sqrt(2.0))*(2*x.^2-1).*exp((-1.0/2.0)*x.^2);

hold on
plot(x,psi_0, '--')
plot(x(1:10:1000),psi_1(1:10:1000), '.-')
plot(x,psi_2)
legend('psi_0','psi_1','psi_2')
xlabel('x')
ylabel('psi_n')
hold off