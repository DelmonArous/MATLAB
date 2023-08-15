x0 = -100.0;                % [fm]
a = 0.0004;                 % [fm^-2]
A = ((2*a)/pi)^(1/4);       % [fm^(-1/2)]

x = linspace(-200,0,1000);  % [fm]

psi0 = A*exp(-a*(x-x0).^2); % [1/fm]

plot(x, abs(psi0).^2)
xlabel('x [fm]')
ylabel('|Psi(x,0)|^2 [1/fm]')