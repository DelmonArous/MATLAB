function SF = MLQ(D, R, lambda0, alpha, beta, delta, str)

if strcmp(str, 'lin')
    SF  = exp(lambda0 + alpha.*D + beta.*(D.^2) + delta.*R);
elseif strcmp(str, 'sq')
    SF  = exp(lambda0 + alpha.*D + beta.*(D.^2) + delta.*(D.*R));
else
    SF  = exp(lambda0 + alpha.*D + delta.*R);
end

end