function [x,y] = newtonfler(x1,y1,N)

x = zeros(1,N);
y = zeros(1,N);
x(1) = x1;
y(1) = y1;
u = [x1;y1];

for n = 1:N
    A = [2*x(n) -1;2*x(n) 2*y(n)];
    v = [(x(n))^2-y(n);x(n)^2+y(n)^2-1];
    u = u-A\v;
    x(n+1) = u(1);
    y(n+1) = u(2);
end
end