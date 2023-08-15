antsec = 3;
fs = 44100;
t = linspace(0, antsec, fs*antsec);

a = 1; 
b = 0.2;

S = a*sin(2*pi*440.*t) + b*sin(2*pi*4400.*t);
S = S/max(abs(S));

plot(t, S)
axis([0 0.01 -1.1 1.1])

% playerobj = audioplayer(S, fs);
% playblocking(playerobj);