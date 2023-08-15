for n = [5 30]
    
    figure
    [X,Y,psi] = streamfun(n);
    contour(X,Y,psi);
    xlabel('x-akse');
    ylabel('y-aske');

    title('Plot av konturlinjer for psi')
    legend(['n = ',num2str(n)])
    
end