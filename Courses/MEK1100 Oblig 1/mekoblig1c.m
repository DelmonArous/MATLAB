t_star = linspace(0,1,1001);
x_star = t_star;
theta = [pi/6 pi/4 pi/3];

for i = 1:3
    if i == 1
        line = '-';
    elseif i == 2
        line = '--';
    else
        line = ':';
    end
    y_star = tan(theta(i))*(-t_star.^2 + t_star);
    hold on
    plot(x_star, y_star,line)
    xlabel('x^*')
    ylabel('y^*')
    hold off
end
legend(['theta = pi/6'],['theta = pi/4'],['theta = pi/3'])