function [x] = newton(x1,N)

x = zeros(1,N);
x(1) = x1;

for n = 1:N
    x(n+1) = x(n) - ((3*x(n)^2 - 1)/(x(n)^3 - x(n)));
    x(n) = x(n+1);
end

disp(x(n));

end