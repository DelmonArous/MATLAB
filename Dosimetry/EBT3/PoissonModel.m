function [lambda] = PoissonModel(D, R, intercept, alpha, beta, eta)

if exist('eta', 'var') || ~isempty(eta)
    % Return modified LQ model
    lambda = exp(intercept + alpha.*D + beta.*(D.^2) + eta.*R);
else
    % Return conventional LQ model
    lambda = exp(intercept + alpha.*D + beta.*(D.^2));
end

end