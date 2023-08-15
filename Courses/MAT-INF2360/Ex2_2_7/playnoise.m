function playnoise(c)

[S fs] = wavread('castanets.wav');

for i = 1:2
   x = S(:,i);
   x = x./max(abs(x));
   N = length(x);
   y = x + transpose(c*(2*rand(1,N) - 1));
   y = y/max(abs(y));
   S(:,i) = y;
end

playerobj = audioplayer(S, fs);
playblocking(playerobj);

end