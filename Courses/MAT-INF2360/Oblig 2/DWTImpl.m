function xnew=DWTImpl(h0,h1,x,m)

N1 = length(h0);
topad1 = (N1-1)/2;
N2 = length(h1);
topad2 = (N2-1)/2;

len = length(x);
for mres=1:m
    xnew = [x((topad1+1):(-1):2); x(1:len); x((len-1):(-1):(len-topad1))];
    x1 = conv(h0,xnew);
    x1 = x1(N1:(length(x1)-(N1-1)));

    xnew = [x((topad2+1):(-1):2); x(1:len); x((len-1):(-1):(len-topad2))];
    x2 = conv(h1,xnew);
    x2 = x2(N2:(length(x2)-(N2-1)));

    % Reorganize the coefficients
    l = x1(1:2:length(x1));
    h = x2(2:2:length(x2));
    x(1:len) = [l; h];
    len = ceil(len/2);
end

xnew = x;

end
  