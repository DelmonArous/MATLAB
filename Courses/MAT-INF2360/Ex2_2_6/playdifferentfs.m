function playdifferentfs()

[S fs] = wavread('castanets.wav');

x = S(:,2);
x = x./max(abs(x));

for i = [fs/2, fs, 2*fs]
    playerobj = audioplayer(S,i);
    playblocking(playerobj);
end

end