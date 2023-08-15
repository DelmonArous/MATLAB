function playDWTfilterslowerdifference(m,h0,h1,g0,g1)

[S fs] = wavread('castanets.wav');
x = S(:,1); % velger en kanal
x = x./max(abs(x)); % verdier mellom -1 og 1
x = x(1:2^17); % velger de 2^17 første lydsamplene av wav-filen

% m'te grads DWT ved å anvende filtrene h0 og h1 m ganger
xnew = DWTImpl(h0,h1,x,m); 
N = length(xnew);

% beholder kun wavelet koeffisienter som representerer detalj
% og setter bidrag fra V0 til null
xnew(1:N/2^m) = 0.;

% IDWT på de resulterende koeffisientene ved å anvende filtrene 
% g0 og g1 m ganger
xnew = IDWTImpl(g0,g1,xnew,m);

% rekonsturerer og spiller av resulterende lyd
playerobj = audioplayer(xnew, fs);
playblocking(playerobj);

end