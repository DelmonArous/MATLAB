x0 = [1; 0; 0];
A = [8 0 12; 1 -2 1; 0 3 0];
n = 5;

for k = 1:n
    x_k = (A^k)*x0;
    my_k = max(x_k);
    x_k1 = (1./my_k)*x_k;
   
    disp('x_k: ')
    disp(x_k)
    disp('my_k: ')
    disp(my_k)
end