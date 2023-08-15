function SF = LQ(D, lambda0, alpha, beta)

SF = exp(lambda0 + alpha.*D + beta.*(D.^2));

end