clc;
clear all;

x = [0.4992 -0.8661 0.7916 0.9107 0.5357...
    0.6574 0.6353 0.0342 0.4988 -0.4607];
alpha0 = 0;

f = @(alpha)(-sum( log( (1+alpha*x)./2 ) ));
df = @(alpha)(-sum( x./(1+alpha*x) ));
d2f = @(alpha)(sum( x.^2./(1+alpha*x).^2 ));

[alphaopt, numit] = newtonbacktrack(f, df, d2f, alpha0);