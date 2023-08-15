function reducebass(k)

[S fs] = wavread('castanets.wav');

x = S(:,1);
x = x./max(abs(x));
N = length(x);

c = [1];
kernel = [1 1];
if (k > 0) 
    for i = 1:2*k
        c = conv(c, kernel);
    end
end

y = zeros(1,N);

for t = (k+1):(N-k)
   for j = 1:(2*k + 1)
       y(t) = y(t) + ((-1)^(k+1-j)*c(j)*x(t+j-k-1))/2^k;
   end
end

y = y./max(abs(y));

playerobj = audioplayer(y, fs);
playblocking(playerobj);

end