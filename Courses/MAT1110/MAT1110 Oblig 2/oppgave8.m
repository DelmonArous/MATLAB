A = zeros(121,121);
b = zeros(121, 1);

for i = 0:10

     for j = 0:10
        k = (i + 1) + 11*j; 
        
        if i == 0
           b(k,1) = j;
        elseif j == 0;
            b(k,1) = i;
        elseif i == 10;
            b(k,1) = 10-j;
        elseif j == 10;
            b(k,1) = 10-i;
        else
            b(k,1) = 0;
        end
        
        
        A(k,k) = 1;
        
        if (0<i && i<10) && (0<j && j<10)
            A(k,k) = -4;
            A(k,k-11) = 1;
            A(k,k-1) = 1;
            A(k,k+1) = 1;
            A(k,k+11) = 1;
        end
        
    end
end

sum(sum(A));
sum(b);
rref([A b]);
x = A\b;
     
[u v] = meshgrid(0:10,0:10);
k = (u+1) + 11*v;
f = x(k);
surf(u,v,f)