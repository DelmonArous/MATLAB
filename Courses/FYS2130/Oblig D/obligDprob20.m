N = 2^16;
Nstart = 1.0;
tmin = 0.0;
tmax = 2.0;
t = linspace(0, tmax*(N-1)/N, N);

s = 'enghornH.wav';
[f, Fs, type] = wavread(s, [Nstart Nstart+N-1]);
g = f(:,1);

T = tmax - tmin;

plot(t, g)
title('Opprinnelig signal')
xlabel('Tid [s]')
ylabel('Relativ amplitude')

z = (1/sqrt(N))*fft(g);
deltaF = 1.0/T;

for i = 1:N
    F(i) = deltaF*(i-1);
end

figure;
plot(F, abs(z), '.-k')
xlabel('Frekvens [Hz]')
ylabel('|X_k|')
