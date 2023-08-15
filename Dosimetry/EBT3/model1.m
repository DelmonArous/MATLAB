function [D] = model1(netOD, a, b, n)

D = (a .* netOD) + (b .* (netOD.^n));

end