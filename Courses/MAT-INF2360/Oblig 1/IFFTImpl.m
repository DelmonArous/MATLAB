% Compute the IDFT of the column vector x using the IFFT algorithm
function y = IFFTImpl(x)

N = length(x);

if (N == 1) 
    y = x(1);
else
    ye = IFFTImpl(x(1:2:(N-1)));
    yo = IFFTImpl(x(2:2:N));
    D = exp(2*pi*1i*(0:(N/2-1))/N);
    yo = D*yo;
    y = [ye + yo; ye - yo]/sqrt(2);
end

end