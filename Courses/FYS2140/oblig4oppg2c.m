A = 1.0;
lambda = 1.0;
sigma = 1/(sqrt(2)*lambda);
x = linspace(-10,10,1000);

psi = A*exp(-lambda*abs(x));
psi_sigma = A*exp(-lambda*abs(sigma));

plot(x,(abs(psi)).^2)
xlabel('x')
ylabel('|psi|^2')
legend(['A=' num2str(A)...
    ', lambda=' num2str(lambda)])
    text(sigma, (abs(psi_sigma))^2,'\leftarrow x1',...
    'HorizontalAlignment','left')
text(-sigma, (abs(psi_sigma))^2,'x2 \rightarrow',...
    'HorizontalAlignment','right')
