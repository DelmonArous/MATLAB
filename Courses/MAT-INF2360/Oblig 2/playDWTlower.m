function playDWTlower(m)

[S fs] = wavread('castanets.wav');
x = S(:,1); % velger en kanal
x = x./max(abs(x)); % verdier mellom -1 og 1
x = x(1:2^17); % velger de 2^17 første lydsamplene av wav-filen

% m'te grads DWT ved bruk av Haar-waveleten med tilhørende filtre h0 og h1
xnew = DWTHaarImpl(x, m);
N = length(xnew);

% setter wavelet koeffisienter som representerer detalj til null
% og beholder kun bidrag fra V0
xnew(1+2^(17-m):N) = 0.;

xnew = IDWTHaarImpl(xnew, m); % IDWT på de resulterende koeffisientene

% rekonsturerer og spiller av resulterende lyd
playerobj = audioplayer(xnew, fs);
playblocking(playerobj);

end