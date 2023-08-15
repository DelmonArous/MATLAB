x0 = [1; 0; 0];
A = [8 0 12; 1 -2 1; 0 3 0];
alpha = -1.4;
n = 5;

for k = 0:n
    y_k = (A-alpha*eye(3))\x0;
    my_k = abs(max(y_k));
    ny_k = alpha + (1./my_k);
    x_k1 = (1/my_k)*y_k;
    
    disp('ny_k: ')
    disp(ny_k)
    disp('x_k: ')
    disp(x_k)
end