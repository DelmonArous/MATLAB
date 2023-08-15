x = linspace(-200, 200, 1000);

for i = 0:10:180
    t = i*10^(-4);
    psi = Psixt(x, t);
    plot(x, psi);
    xlabel('x [fm]')
    ylabel('|Psi(x,t)|^2 [1/fm]')
    %axis([-200 200 0 0.016])
    legend(['t=',num2str(t),' as'])
    drawnow;
    pause(0.005);
end 