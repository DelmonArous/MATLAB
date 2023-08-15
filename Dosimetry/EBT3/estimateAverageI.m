function [I_avg, sigma_avg] = estimateAverageI(I, sigma)

N = length(I);
w           = (1./(sigma.^2)) ./ (sum(1./(sigma.^2)));
I_avg       = sum(w.*I);
sigma_avg   = sqrt(N ./ (sum(1./(sigma.^2))));

end