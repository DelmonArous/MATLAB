x = linspace(-8, 8, 1000);
l = linspace(0.445,0.447,21);
sum = zeros(length(l),1);         

for j= 1:length(l)
    sumprob = 0;
    for i = 0:10:200;
        t = i*10^-4;
        psixtsq = Psixt_v2(x,t,l(j)); % definerer en ny funksjon 
                                      % som tar inn tre argumenter
        prob = trapz(x, psixtsq);
        sumprob = sumprob + prob;
    end
    % Sannsynlighet for at partikkelen befinner seg 
    % i bronnen for l(j)
    sum(j)=sumprob;
end
 maxprob = max(sum);    % finner maksimal sannsynlighet for at
                        % partikkelen befinner seg i bronnen            
% Finner l som maksimerer sannsynligheten for at
% partikkelen befinner seg i bronnen for alle tidene
for i = 2:length(sum)
    if sum(i-1)<sum(i)
        l(i);
    end
end