function [sigma] = sigma_tot(sigma_mono, sigma_E, alpha, p, E)

sigma = sqrt(sigma_mono^2 + sigma_E^2*alpha^2*p^2*E^(2*p-2));

end