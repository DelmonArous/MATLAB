function [u] = PoissonDirichletSolver(f, n)
    h = 1./(n+1);
    u = zeros(1, n+2);
    alpha = zeros(1, n+2);
    beta = zeros(1, n+2);
    
    for i = 1:(n+1)
        x_i = (i-1)*h;
        alpha(i+1) = alpha(i) + (h/4.)*( f(x_i) + ...
            2*f(x_i + h/2.) + f(x_i+h) );
        beta(i+1) = beta(i) + (h/4.)*(x_i*f(x_i) + ...
            2.*(x_i+h/2.)*f(x_i + h/2.) + (x_i+h)*f(x_i +h) );
    end
    
    for i = 2:(n+1)
        x_i = (i-1)*h;
        u(i) = x_i*(alpha(n+2) - beta(n+2)) + beta(i) - x_i*alpha(i);
    end
end