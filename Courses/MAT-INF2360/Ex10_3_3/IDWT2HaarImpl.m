function X = IDWT2HaarImpl(Y, m)

for mres = 1:m
    l1 = size(Y,1)/2^(mres-1);
    l2 = size(Y,2)/2^(mres-1);
    
    for s = 1:l2
        Y(1:l1,s) = IDWTHaarImpl(Y(1:l1,s),1);
    end
    
    Y = Y';
    
    for s = 1:l1
        Y(1:l2,s) = IDWTHaarImpl(Y(1:l2,s),1);
    end
    
    Y = Y';
end

X = Y;

end