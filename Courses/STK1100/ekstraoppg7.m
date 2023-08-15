clear all;

gutter = 0:12;
jenter = 12:-1:0;
familier = [3 24 104 286 670 1033 1343 1112 829 478 181 45 7];

n = 12;

p1 = 0.50;
p_k1 = [];
p2 = 0.52;
p_k2 = [];

rel_frek = [];

for k = 1:13
    rel_frek(k) = familier(k)/sum(familier);
    p_k1(k) = binopdf(k-1,n,p1);
    p_k2(k) = binopdf(k-1,n,p2);
end

figure
plot(gutter,p_k1, 'b')
hold on
plot(gutter,p_k2, 'r')
plot(gutter, rel_frek, '--b')
ylabel('Sannsynlighet')
xlabel('Antall gutter')
legend('p=0.50','p=0.52','relativ frekvens')