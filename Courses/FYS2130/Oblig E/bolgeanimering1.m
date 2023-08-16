clear all;

A = 1.5;
k = 0.5;
w = 16.0;
N = 1000;
x = linspace(0,30,N);
y = linspace(0,40,N);
p = plot(x,y,'-','EraseMode','xor');
xlabel('Posisjon [m]')
ylabel('Relativ amplitude [m]')
legend(['A=',num2str(A),', k=',num2str(k),...
    ', w=',num2str(w)])

for i = 1:400
    t = i*0.01;
    %y = A*sin(k*x-w*t);
    y = A*sin(k*x-w*t) + A*sin(k*x+0.1*w*t);
    set(p,'Xdata',x,'YData',y)
    drawnow
    pause(0.02);
end