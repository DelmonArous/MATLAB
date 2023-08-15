function playDWTlower(m)

[S fs] = wavread('castanets.wav');
x = S(:,1);
x = x./max(abs(x));
x = x(1:2^17);

xnew = DWTHaarImpl(x, m);
N = length(xnew);

xnew(1+2^(17-m):N) = 0.;

xnew = IDWTHaarImpl(xnew, m);

playerobj = audioplayer(xnew, fs);
playblocking(playerobj);

end