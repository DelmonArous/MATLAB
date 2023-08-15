clear all;
clc;

m = 2;
g0 = [1/2 1 1/2]/sqrt(2);
g1 = [1]/sqrt(2);

xnew = IDWTImpl(g0, g1,...
    [-alpha; -beta; -delta; 0; 0; 0; 0; -gamma;...
    1; 0; 0; 0; 0; 0; 0; 0], m);

N = length(xnew);

g1 = [xnew((N-2):N) xnew(1:N-3)]; % hvorfor blir xnew symmetrisk?
omega = linspace(0,2*pi,1000);
plot(omega,g1(5)+g1(6)*2*cos(omega)+g1(7)*2*cos(2*omega)...
    +g1(8)*2*cos(3*omega)+g1(9)*2*cos(4*omega))

alpha=1/(g0(2)*g1(5)+2*g0(3)*g1(6));
h0=alpha*(-1).^(1:(length(g1))).*g1;
h1=alpha*(-1).^(0:(length(g0)-1)).*g0;

for m = 1:4
    playDWTfilterslower(m,h0,h1,g0,g1);
    playDWTfilterslowerdifference(m,h0,h1,g0,g1);
end
