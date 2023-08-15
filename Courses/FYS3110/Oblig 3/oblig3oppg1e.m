x = linspace(-10,10,100);
P_B = (sign(x) - sqrt(1 + 8*x.^2)).^2./(8*x.^2 + (sign(x) - sqrt(1 + 8*x.^2)).^2)

plot(x, P_B)
title('The probability P_B as a function of x=g/E_B')
xlabel('x=g/E_B')
ylabel('P_B')