function [x,y] = oppgave551(x1,y1,N)

x = zeros(1,N);
y = zeros(1,N);
x(1) = x1;
y(1) = y1;

for n = 1:N-1
    x(n+1) = 0.5*sin(x(n)+y(n));
    y(n+1) = 0.5*cos(x(n)+y(n));
end
    
end