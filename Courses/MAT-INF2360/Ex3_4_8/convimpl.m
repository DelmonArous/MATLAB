function z = convimpl(x,y)

N = length(x);
M = length(y);

x = [x, zeros(1,M)];
y = [y, zeros(1,N)];

z = zeros(1,M+N-1);

for n = 1:length(z)
    for k = 1:N
        if (n-k+1 > 0)
            z(n) = z(n) + x(k)*y(n-k+1);
        end
    end
end

end