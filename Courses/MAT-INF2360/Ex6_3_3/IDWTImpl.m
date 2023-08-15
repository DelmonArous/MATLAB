function x=IDWTImpl(g0,g1,xnew,m)

  % function changecolumnrows begin
  % Keep track of the indices of the biggest even and odd index.
  g0=g0((length(g0)+1)/2:length(g0));
  g1=g1((length(g1)+1)/2:length(g1));
  maxeveng0=length(g0);
  maxoddg0=length(g0);
  if mod(maxeveng0,2)==0 
      maxeveng0=maxeveng0-1; 
  else
      maxoddg0=maxoddg0-1;
  end
  
  maxeveng1=length(g1);
  maxoddg1=length(g1);
  if mod(maxeveng1,2)==0 
      maxeveng1=maxeveng1-1; 
  else
      maxoddg1=maxoddg1-1;
  end;
  
  a0length=max(maxeveng0,maxoddg1);
  a1length=max(maxeveng1,maxoddg0);
  
  a0=zeros(1,a0length);
  a0(1:2:maxeveng0)=g0(1:2:maxeveng0);
  a0(2:2:maxoddg1)=g1(2:2:maxoddg1);
  
  a1=zeros(1,a1length);
  a1(1:2:maxeveng1)=g1(1:2:maxeveng1);
  a1(2:2:maxoddg0)=g0(2:2:maxoddg0);
  
  a0=[a0(length(a0):(-1):2) a0];
  a1=[a1(length(a1):(-1):2) a1];
  % function changecolumnrows end
  
  N1=length(a0);
  topad1=(N1-1)/2;
  N2=length(a1);
  topad2=(N2-1)/2;
  
  lenall=zeros(m+1,1);
  lenall(1)=length(xnew);
  for mres=1:m
      lenall(mres+1)=ceil(lenall(mres)/2);
  end
  % lenall is now the lengths of the (lowpass,highpass) sections of xnew
  
  for mres=m:(-1):1
    
    % Reorganize the coefficients first
    l=xnew(1:lenall(mres+1));
    h=xnew((lenall(mres+1)+1):lenall(mres));
    
    xnew(1:2:lenall(mres))=l;
    xnew(2:2:lenall(mres))=h;
    
    x=[xnew((topad1+1):(-1):2); xnew(1:lenall(mres)); xnew((lenall(mres)-1):(-1):(lenall(mres)-topad1))];
    x1=conv(a0,x);
    x1=x1(N1:(length(x1)-(N1-1)));
  
    x=[xnew((topad2+1):(-1):2); xnew(1:lenall(mres)); xnew((lenall(mres)-1):(-1):(lenall(mres)-topad2))];
    x2=conv(a1,x);
    x2=x2(N2:(length(x2)-(N2-1)));
    
    xnew(1:2:length(x1))=x1(1:2:length(x1));
    xnew(2:2:length(x2))=x2(2:2:length(x2));
  end
  x=xnew;
end