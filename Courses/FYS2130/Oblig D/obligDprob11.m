clear all;

tmin = 0.0;
tmax = 15.0;
N = 512;
t = linspace(tmin, tmax*(N-1)/N, N);

f = sin(((2*pi*13.2)/tmax)*t);
%f = 0.7*sin(2*pi*50*t) + cos(2*pi*120*t);

plot(t, f)
title('Opprinnelig signal')
xlabel('Tid [s]')
ylabel('Relativ amplitude')

z = (1/sqrt(N))*fft(f);
deltaF = 1.0/(tmax-tmin);

for i = 1:N
    F(i) = deltaF*(i-1);
end

figure;
plot(F, abs(z), '.-k')
xlabel('Frekvens [Hz]')
ylabel('|X_k|')