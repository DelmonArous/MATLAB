m = 100;                % antall kast

Y = unidrnd(2, [1, m]);
X = (Y == 2);           % true for Y == 2, false ellers

k= sum(X);              % antall kron i m kast
p = k/m;                % gir 0.5

K = cumsum(X);
M = 1:m;
Q = K./M;

for i = [10,30,100]
    q = sum(X(1:i))/i;  % gir normalfordelt om 0.5
    K(i) == sum(X(1:i));
    K(i) == sum(X(1:i-1)) + X(i);
    K(i) == K(i-1) + X(i);
end

plot(M,Q)
title('Anslag myntkast Q som funksjon av antall myntkast M');
xlabel('M');
ylabel('Q');