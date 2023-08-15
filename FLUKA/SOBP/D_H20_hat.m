function [D] = D_H20_hat(Phi0, R0, epsilon, z)

D = (Phi0/(1+0.012*R0)) .* (17.93.*(R0-z).^(-0.435) + ...
    (0.444+31.7*epsilon/R0).*(R0-z)^(0.565) );

end