function [psixt] = Psixt(x,t)
% Funksjon som evaluerer |Psi(x,t)|^2 for endelig-brønn potensialet

% Definer konstanter
% Enheter er fm (lengde), MeV (energi) og as (tid)
a = .0004; mc2 = 939.6; hbarc = 197.3; c = 3.0*10^5;
x0 = -100.0; l = 0.2; R = 8.0; V0 = 35.0;

% Integral over bølgetall
% Tabell med k-verdier
ki = l-0.3; kf = l+0.3; nk = 1000;
k = ki : (kf-ki)/nk : kf;
% På grunn av faktoren exp(-(k(ni)-l)^2/4/a) i \phi kan vi begrense
% integrasjonsområdet til k ganske sterkt uten å miste noe presisjon

% Integrasjon over k
sum = 0;
for ni = 1 : length(k)     % Loop over k
  % Vi bruker her den \phi(k) som vi fant et eksplisitt uttrykk for for
  % den frie partikkelen. Hvorfor kan vi gjøre det? Vi burde jo egentlig
  % finne en ny \phi(k) uttrykt ved hjelp av løsningene for den endelige
  % brønnen. Saken er at bølgefunksjonen begynner sentrert om et området 
  % der partikkelen er fri, og et stykke fra effekten av brønnen slik at 
  % dette er en god approksimasjon. Hvor god kan testes ved å se på plot(x,Psixt(x,0));
  % Konsekvensen av dette valget er noe numerisk ustabilitet i rutinen, 
  % spesielt hvis de valgte parameterene, slik som l, endres dramatisk.
  % Kanskje du kan forbedre denne rutinen ved å finne den riktige \phi(k) 
  % numerisk eller analytisk?
  phik = 1/(2*pi*a)^(1/4)*exp(-1i*x0*(k(ni)-l))*exp(-(k(ni)-l)^2/4/a);
  omega = hbarc*k(ni)^2/2/mc2*c;
  
  % Her må vi sette inn løsninger for den endelige brønnen
  kappa0 = 2*mc2*V0/hbarc^2;
  kappa = sqrt(k(ni)^2+kappa0);
  A = 1;  % A setter normaliseringen på bølgefunksjonen, som vi ignorerer her
  F = exp(-2i*k(ni)*R)/(cos(2*kappa*R)-1i*(k(ni)^2+kappa^2)/k(ni)/kappa*sin(2*kappa*R))*A;
  B = 1i*sin(2*kappa*R)/(2*k(ni)*kappa)*kappa0*F;
  C = (sin(kappa*R)+1i*k(ni)/kappa*cos(kappa*R))*exp(1i*k(ni)*R)*F;
  D = (cos(kappa*R)-1i*k(ni)/kappa*sin(kappa*R))*exp(1i*k(ni)*R)*F;
  
  psix = (A*exp(1i*k(ni)*x)+B*exp(-1i*k(ni)*x)).*(x < -R) + ...
   + (C*sin(kappa*x)+D*cos(kappa*x)).*(-R <= x & x <= R) + ...
   + (F*exp(1i*k(ni)*x)).*(x > R);
  
  sum = sum + phik*psix*exp(-1i*omega*t);
end % Slutt k-loop

psixt = sum*(kf-ki)/nk;
psixt = psixt/sqrt(2*pi);
psixt = abs(psixt).^2;

end