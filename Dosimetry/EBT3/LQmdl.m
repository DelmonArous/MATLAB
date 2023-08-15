function [mu, sigma] = LQmdl(D, lambda0, alpha, beta, ...
    sigma_lambda0, sigma_alpha, sigma_beta)

mu    = exp(lambda0 + alpha.*D + beta.*(D.^2));
sigma = sqrt( (mu.*sigma_lambda0).^2 + (D.*mu.*sigma_alpha).^2 + ...
    ((D.^2).*mu.*sigma_beta).^2 ); 

end