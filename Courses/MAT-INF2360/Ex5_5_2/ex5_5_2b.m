clear all;
clc;

A = zeros(4,4);
e = zeros(4,1);

for k = 0:3
    A(k+1,:) = [quad(@(t)t.^k.*(1-abs(t)),-1,1) ...
        quad(@(t)t.^k.*(1-abs(t-1)),0,2) ...
        quad(@(t)t.^k.*(1-abs(t+1)),-2,0) ...
        quad(@(t)t.^k.*(1-abs(t-2)),1,3)];
    e(k+1) = quad(@(t)t.^k.*(1-2*abs(t-0.5)),0,1);
end

x = A\e;
alpha = x(1);
beta = x(2);
gamma = x(3);
delta = x(4);

t = linspace(-2, 3, 1000);
psi = (1/sqrt(2))*(t>=0).*(t<=1).*(1-2*abs(t-0.5)) - ...
    alpha*(t>=-1).*(t<=1).*(1-abs(t)) - ...
    beta*(t>=0).*(t<=2).*(1-abs(t-1)) - ...
    gamma*(t>=-2).*(t<=0).*(1-abs(t+1)) - ...
    delta*(t>=1).*(t<=3).*(1-abs(t-2));

plot(t, psi)
