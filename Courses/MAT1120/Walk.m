function [x2,x3,x4] = Walk(p2,p3,p4)

if (0<p2 && p2<1 && 0<p3 && p3<1 && 0<p4 && p4<1)
    q2 = 1 - p2;
    q3 = 1 - p3;
    q4 = 1 - p4;
    
    A = [1 -q2 0; -p3 1 -q3;...
         0 -p4 1];
    b = [p2; 0; 0];
    y = A\b;
    x2 = y(1);
    x3 = y(2);
    x4 = y(3);
else
    disp('Must have 0 < p_j < 1');
end
end

% Calling on the function with given values
% (p2,p3,p4) = (1/2,2/3,1/4) gives
% x2 = 0.7857 
% x3 = 0.5714
% x4 = 0.1429