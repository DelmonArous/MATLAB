function y=DFTImpl(x)

N = length(x);
FN = zeros(N);

for n = 1:N
    FN(n,:) = exp(-2*pi*1i*(n-1)*(0:(N-1))/N)/sqrt(N);
end

y = FN*x;

end