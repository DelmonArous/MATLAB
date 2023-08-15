clear all;
n = [3 10 30];
mu = [0.5 1.0 0.5];
sigma = [1/(2*sqrt(3)) 1.0 0.5];
int = [-Inf, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, Inf];
Z_normcdf = [0.0062 0.0165 0.0441 0.0918 0.1499 0.1915...
    0.1915 0.1499 0.0918 0.0441 0.0165 0.0062];

% for i = 1:3
%     X_unifrnd = unifrnd(0,1,n(i),10000);
%     meanX= mean(X_unifrnd);
%     Z = sqrt(n(i))*(meanX-mu(1))/sigma(1);   
%     figure
%     hist(Z,-3:0.25:3)
%     xlabel('z_n')
%     ylabel('Hyppighet')
%     title({'Ikke-normert histogram av av de standardiserte gjennomsnittene Z_n';...
%         ['for den uniforme fordelingen; n=', num2str(n(i))]})
%     
%     figure
%     hold on
%     ant = histc(Z, int);
%     relfrekv = ant(1:12)/10000;
%     length(relfrekv)
%     plot(-5:6,relfrekv), plot(-5:6,Z_normcdf,'--')
%     legend('Relativ frekvens','\Phi(z)')
%     xlabel('z')
%     ylabel('Sannsynlighet')
%     title(['Plot av relativ frekvens for Z_n og \Phi(z); n=', num2str(n(i))])
%     hold off
% end  

for i = 1:3
    X_exprnd = exprnd(1,n(i),10000);
    meanX= mean(X_exprnd);
    Z = sqrt(n(i))*(meanX-mu(2))/sigma(2);
    figure
    hist(Z,-3:0.25:3)
    xlabel('z_n')
    ylabel('Hyppighet')
    title({'Ikke-normert histogram av av de standardiserte gjennomsnittene Z_n';...
        ['for eksponensialfordelingen; n=', num2str(n(i))]});
    
    figure
    hold on
    ant = histc(Z, int);
    relfrekv = ant(1:12)/10000;
    plot(-5:6,relfrekv), plot(-5:6,Z_normcdf,'--')
    legend('Relativ frekvens','\Phi(z)')
    xlabel('z')
    ylabel('Sannsynlighet')
    title(['Plot av relativ frekvens for Z_n og \Phi(z); n=', num2str(n(i))])
    hold off
end 

for i = 1:3
    X_binornd = binornd(1,0.5,n(i),10000);
    meanX= mean(X_binornd);
    Z = sqrt(n(i))*(meanX-mu(3))/sigma(3);
    figure
    hist(Z,-3:0.25:3)
    xlabel('z_n')
    ylabel('Hyppighet')
    title({'Ikke-normert histogram av av de standardiserte gjennomsnittene Z_n';...
        ['for Bernoulli-fordelingen; n=', num2str(n(i))]})

    figure
    hold on
    ant = histc(Z, int);
    relfrekv = ant(1:12)/10000;
    plot(-5:6,relfrekv), plot(-5:6,Z_normcdf,'--')
    legend('Relativ frekvens','\Phi(z)')
    xlabel('z')
    ylabel('Sannsynlighet')
    title(['Plot av relativ frekvens for Z_n og \Phi(z); n=', num2str(n(i))])
    hold off
end 