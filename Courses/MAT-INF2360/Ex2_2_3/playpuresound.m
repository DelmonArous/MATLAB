function playpuresound(f)

antsec = 3;
fs = 2.5*f;
t = linspace(0, antsec, fs*antsec);
S = sin(2*pi*f*t);

playerobj = audioplayer(S, fs);
playblocking(playerobj);
wavwrite(S, fs, 'puretonef.wav');

end