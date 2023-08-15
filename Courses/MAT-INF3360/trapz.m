function [I] = trapz(F, a, b, n)
    h = (b-a)/(n+1);
    I = F(a)/2. + F(b)/2.;
    
    for i = 1:n
        x = i*h;
        I = I + F(x);
    end
    
    I = I*h;
end