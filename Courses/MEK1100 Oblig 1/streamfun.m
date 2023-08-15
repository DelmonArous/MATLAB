% Regner ut et grid og en stromfunksjon
function [X,Y,psi] = streamfun(n)

%setter standardverdi for n til 20
if nargin < 1; 
    n=20;
end

%linspace brukes istedet for  -0.5*pi:step:0.5*pi for å unngå
%toys pga. avrunding
x=linspace(-0.5*pi,0.5*pi,n);
%resultatet er en vektor med n elementer, fra -pi/2 til pi/2
[X,Y] = meshgrid(x,x);
psi=cos(X).*cos(Y);

end