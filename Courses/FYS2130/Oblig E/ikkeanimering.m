clear all;

A = 1.5;
k = 0.5;
w = 16.0;
N = 1000;
x = linspace(0,30,N);
t = 50.0;

%y = A*sin(k*x-w*t);
y = A*sin(k*x-w*t) + A*sin(k*x+w*t);

plot(x,y,'-','EraseMode','xor');
xlabel('Posisjon [m]')
ylabel('Relativ amplitude [m]')
legend(['A=',num2str(A),', k=',num2str(k),...
    ', w=',num2str(w)])
