function lambda = lambda_n(n)

lambda = (1 + cos(pi*n/2) + cos(pi*n/4)...
    + cos(pi*n/4)*cos(pi*n/2));

end