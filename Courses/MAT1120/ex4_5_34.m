C = [1 0 -1 0 1 0 -1; 0 1 0 -3 0 5 0;...
    0 0 2 0 -8 0 18; 0 0 0 4 0 -20 0;...
    0 0 0 0 8 0 -48; 0 0 0 0 0 16 0;...
    0 0 0 0 0 0 32];

rref(C)

for i = 1:7
    A(:,i) = cos(t).^i;
end