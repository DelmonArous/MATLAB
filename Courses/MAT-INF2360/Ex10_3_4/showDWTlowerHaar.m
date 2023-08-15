function showDWTlowerHaar(m)

X = double(imread('lena.png','png'));

for k = 1:3
   X(:,:,k) = DWT2HaarImpl(X(:,:,k),m); 
end

[l1,l2,l3] = size(X);
tokeep = X(1:(l1/2^m),1:(l2/2^m),:);
X = zeros(size(X));
X(1:(l1/2^m),1:(l2/2^m),:) = tokeep;

for k = 1:3
    X(:,:,k) = IDWT2HaarImpl(X(:,:,k),m);
end

imshow(uint8(255*mapto01(X)));

end