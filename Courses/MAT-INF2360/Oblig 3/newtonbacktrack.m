function [xopt,numit] = newtonbacktrack(f,df,d2f,x0)
  beta = 0.2; s = 0.5; sigma = 10^(-3);
  epsilon = 10^(-3);
  xopt = x0;
  maxit = 100;
  
  for numit = 1:maxit
      d = -d2f(xopt)\df(xopt);
      eta = -df(xopt)'*d;
      if (eta/2 < epsilon)
          break;
      end
      m = 0;
      while (f(xopt)-f(xopt+beta^m*s*d) < -sigma *beta^m*s *(df(xopt))'*d)
        m = m+1;
      end
      xopt = xopt+beta^m*s*d;
  end
  
end