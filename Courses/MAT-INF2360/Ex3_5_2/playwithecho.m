function playwithecho(c,d)

[S fs] = wavread('castanets.wav');

x = S(:,1);
x = x./max(abs(x));
N = length(x);

y = [transpose(x) zeros(1,d)];
y((d+1):N) = x((d+1):N) + c*x(1:(N-d));

playerobj = audioplayer(y, fs);
playblocking(playerobj);

end