V0 = 35.0;              % [MeV]
R = 8.0;                % [fm]
m = 939.6;              %/c^2 [MeV]
hbarc = 1240/(2*pi);    % [fm MeV]

E = linspace(0, 150, 1000); % [MeV] 

T = 1./(1+(V0^2./(4.*E.*(E+V0))).*(sin((2*R/hbarc)*sqrt(2*m.*(E+V0)))).^2);
    % ingen benevning (sannsynlighet for transmisjon)

plot(E, T)
xlabel('E [MeV]')
ylabel('Transmisjonskoeffisient T(E)')
axis([0 150 0 1.1])