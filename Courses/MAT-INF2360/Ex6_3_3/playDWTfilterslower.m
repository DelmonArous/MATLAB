function playDWTfilterslower(m,h0,h1,g0,g1)

[S fs] = wavread('castanets.wav');
x = S(:,1);
x = x./max(abs(x));
x = x(1:2^17);

xnew = DWTImpl(h0,h1,x,m);
N = length(xnew);

xnew(1+2^(17-m):N) = 0.;

xnew = IDWTImpl(g0,g1,xnew,m);

playerobj = audioplayer(xnew, fs);
playblocking(playerobj);

end