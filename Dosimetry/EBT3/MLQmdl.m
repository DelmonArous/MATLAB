function [mu, sigma] = MLQmdl(D, R, lambda0, alpha, beta, delta, ...
    sigma_lambda0, sigma_alpha, sigma_beta, sigma_delta)

mu    = exp(lambda0 + alpha.*D + beta.*(D.^2) + delta.*D.*R);
sigma = sqrt( (mu.*sigma_lambda0).^2 + (D.*mu.*sigma_alpha).^2 + ...
    ((D.^2).*mu.*sigma_beta).^2 + (D.*R.*mu.*sigma_delta).^2 ); 

end