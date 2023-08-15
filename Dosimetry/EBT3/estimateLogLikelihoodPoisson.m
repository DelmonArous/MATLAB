function [loglik] = estimateLogLikelihoodPoisson(y, yhat)

loglik = sum(y.*log(yhat) - yhat - log(factorial(y)));

end