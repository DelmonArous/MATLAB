function [sigma] = sigma_mono(alpha_prime, alpha, p, R)

sigma = sqrt(alpha_prime*((p^3*alpha^(2/p))/(3*p-2))*R^(3-2/p));

end