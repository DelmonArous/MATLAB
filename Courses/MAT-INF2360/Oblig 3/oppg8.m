clc;
clear all;
close all;

m = 1;
alpha0 = 0;
a = 0.5;

lineSpec = {'b-','r-','g-'};
legendName = cell(3,1);

teller = 1;
for n = [10 1000 100000]
    % Nullstiller maximum-likelihood-vektoren 
    % for et nytt sett med maaleverdier
    alphaverdier = zeros(10,1);
    
    % Finner maximum-likelihood-estimatet 10 ganger 
    % for hvert sett med n maaleverdier
    for i = 1:10
        x = randmuon(a, m, n); % genererer n maaleverdier for x
        f = @(alpha)(-sum( log( (1+alpha*x)./2 ) ));
        df = @(alpha)(-sum( x./(1+alpha*x) ));
        d2f = @(alpha)(sum( x.^2./(1+alpha*x).^2 ));
        [alphaopt, numit] = newtonbacktrack(f, df, d2f, alpha0);
        alphaverdier(i) = alphaopt;
    end
    
    hold on
    plot(alphaverdier, lineSpec{teller})
    legendName{teller} = ['n = ', num2str(n)];
    teller = teller + 1;
end
ylabel('\alpha');
legend(legendName);