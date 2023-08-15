function ret = randmuon(alpha,m,n)
  
ret = (-1+sqrt(1-alpha*(2-alpha-4*rand(m,n))))/alpha;
  
end