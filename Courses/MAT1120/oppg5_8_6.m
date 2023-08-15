x0 = [0; 1];
A = [6 5; 1 2];
n = 5;

for k = 1:n
    x_k = (A^k)*x0;
    my_k = max(x_k);
    x_k1 = (1/my_k)*x_k;
    
    disp('x_k: ')
    disp(x_k1)
    disp('my_k: ')
    disp(my_k)
end
