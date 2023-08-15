function [F] = MLQ_twocomp(params, D_p, D_v, G_p, G_v, SF_p, SF_v)

alpha   = params(1);
beta    = params(2);
delta_p = params(3);
delta_v = params(4);
% lambda  = params(5);

f_p = SF_p - exp(-alpha.*D_p - beta.*(D_p.^2)); % + delta_p.*G_p);
f_v = SF_v - exp(-alpha.*D_v - beta.*(D_v.^2)); % + delta_v.*G_v);

F = [f_p, f_v];

end