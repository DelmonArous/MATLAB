% Constants
lambda = 1.0;
a = 1.0;
A = sqrt(lambda/pi);
x = linspace(-10,10,1000);
% Gaussian distribution
rho = A*exp(-lambda*(x-a).^2);
% Plot
plot(x, rho)
xlabel('x')
ylabel('rho')
legend(['lambda=' num2str(lambda) ', a=' num2str(a)])