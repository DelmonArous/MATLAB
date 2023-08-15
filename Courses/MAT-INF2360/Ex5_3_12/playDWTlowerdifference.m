function playDWTlowerdifference(m)

[S fs] = wavread('castanets.wav');
x = S(:,1);
x = x./max(abs(x));
x = x(1:2^17);

xnew = DWTHaarImpl(x, m);
N = length(xnew);

xnew(1:N/2^m) = 0.;

xnew = IDWTHaarImpl(xnew, m);

playerobj = audioplayer(xnew, fs);
playblocking(playerobj);

end