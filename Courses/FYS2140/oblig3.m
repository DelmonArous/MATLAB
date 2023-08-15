k0 = 0.65;
k1 = 0.6;
k2 = 0.7;
dk = 0.05;
x = linspace(-50,50,1001);

omega0 = sqrt((k0^2)+1);
omega1 = sqrt((k1^2)+1);
omega2 = sqrt((k2^2)+1);

for t = [0,100]
    y1 = sin(k1*x-omega1*t);
    y2 = sin(k2*x-omega2*t);
    y = y1 + y2;
    y_approx = 2*sin(k0*x-omega0*t).*cos(dk*x-(k0*dk/omega0)*t);
    
    figure;
    plot(x, y, 'b')
    legend(['t=' num2str(t)])
    xlabel('x')
    ylabel('Psi')
    figure;
    plot(x,y-y_approx, 'b')
    legend(['t=' num2str(t)])
    xlabel('x')
    ylabel('Psi')
end