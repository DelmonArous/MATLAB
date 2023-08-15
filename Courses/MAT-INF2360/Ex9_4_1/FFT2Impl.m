function Xnew = FFT2Impl(X)

for k = 1:2
    for s = 1:size(X,2)
        X(:,s) = FFTImpl(X(:,s));
    end
    X = X';
end

Xnew = X;

end