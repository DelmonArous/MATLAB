clear all;
clc;

A = zeros(2,2);
e = zeros(2,1);

for k = 0:1
    A(k+1,:) = [quad(@(t)t.^k.*(1-abs(t)),-1,1) ...
        quad(@(t)t.^k.*(1-abs(t-1)),0,2)];
    e(k+1) = quad(@(t)t.^k.*(1-2*abs(t-0.5)),0,1);
end

x = A\e;