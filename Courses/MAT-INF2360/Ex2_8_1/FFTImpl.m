function y = FFTImpl(x)

N = length(x);

if N == 1
    y = x(1);
else
    ye = FFTImpl(x(1:2:(N-1)));
    yo = FFTImpl(x(2:2:N));
    D = exp(-2*pi*1j*(0:(N/2-1))/N);
    yo = yo.*D;
    y = [ ye + yo; ye - yo]/sqrt(2);

end