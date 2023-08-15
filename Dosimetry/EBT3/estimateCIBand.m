function [struct] = estimateCIBand(y, x, significance_level)

% Estimate CI range
alpha = [significance_level/2 1-(significance_level/2)]; 

% Estimate CI
N       = size(y, 1);
y_mean  = mean(y);
y_SEM   = std(y) / sqrt(N);

CI95        = tinv(alpha, N-1);
y_CI95      = bsxfun(@times, y_SEM, CI95(:));
CI95_band   = y_mean + y_CI95;

x_CI95      = [x x(end:-1:1)];
CI95_band   = [CI95_band(1,:) fliplr(CI95_band(2,:))];

% Return
struct.y_avg        = y_mean;
struct.x            = x;
struct.CI95_value   = CI95_band;
struct.CI95_x       = x_CI95;

end