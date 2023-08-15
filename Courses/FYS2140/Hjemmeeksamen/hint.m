% Hjemmeeksamen
clear all; 
close all; 
format long;

% Tabell med x-verdier
x0 = -100; xf = 200; xi=-xf;  nx = 1000;
x = xi : (xf-xi)/nx : xf; 
% Tabell med t-verdier
ti = 0; tf = 0.020; nt = 1000;
t = ti : (tf-ti)/nt : tf;

% Definer konstanter
% Enheter er fm (lengde), MeV (energi) og as (tid)
a = .0004; mc2 = 939.6; l = 0.4; hbarc = 197.3; c = 3.0*10^5;

% Bølgefunksjon for t=0
Psi0sq = sqrt(2*a/pi)*exp(-2*a*(x-x0).^2);

% Tegn bølgefunksjon som funksjon av tid
figure;
p = plot(x,Psi0sq);
axis([xi xf 0 0.02]);
xlabel('x [fm]');
ylabel('|Psi(x,t)|² [1/fm]');
T0 = 2*hbarc*a*c/mc2;
fprintf('T_0 = 2*hbar*a/m: %d ',T0); fprintf('\n');
for ni = 1 : length(t)     % Loop over tid
  T = T0*t(ni);
  w = sqrt(a/(1+T^2));
  Psisq = sqrt(2/pi)*w*exp(-2*w^2*(x-x0-T*l/(2*a)).^2);
  set(p,'XData',x,'YData',Psisq)
  drawnow
  pause(0.005);
end % Slutt tidsloop
