g = linspace(0,10,500);

for E_B = [-10 0 10]
    figure;
    lambda_0 = 0.*g;
    lambda_plus = (E_B + sqrt(E_B^2 + 8*g.^2))./2;
    lambda_minus = (E_B - sqrt(E_B^2 + 8*g.^2))./2; 
    hold on
    plot(g, lambda_0, 'b', g, lambda_plus, 'r', g, lambda_minus, 'g')
    legend('\lambda_0', '\lambda_+', '\lambda_-')
    title(['Energy level diagram with eigenvalues as function of g, E_B=', num2str(E_B)])
    xlabel('g');
    ylabel('\lambda_n(g)');
    hold off
end