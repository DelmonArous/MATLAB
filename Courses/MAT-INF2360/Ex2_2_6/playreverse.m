function playreverse()

[S fs] = wavread('castanets.wav');
sz = size(S,1);
newS = [S(sz:(-1):1,1) S(sz:(-1):1,2)];

x = S(:,2);
x = x./max(abs(x));

for i = [fs]
    playerobj = audioplayer(newS,i);
    playblocking(playerobj);
    wavwrite(newS,i,'castanetsreverse.wav');
end
    
end