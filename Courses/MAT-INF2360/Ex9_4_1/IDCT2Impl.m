function X = IDCT2Impl(Xnew)

for k = 1:2
    for s = 1:size(Xnew,2)
        Xnew(:,s) = IDCTImpl(Xnew(:,s));
    end
    Xnew = Xnew';
end

X = Xnew;

end