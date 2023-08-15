[alder menn kvinner] = textread('dod.txt', '', 'headerlines', 7);

qmenn = menn./1000;
qkvinner = kvinner./1000;
pmenn = 1 - qmenn;
pkvinner = 1 - qkvinner;

smenn = cumprod(pmenn);
skvinner = cumprod(pkvinner);
smenn = [1;smenn(1:99)];
skvinner =[1;skvinner(1:99)];

pr_menn = qmenn.*smenn;
pr_kvinner = qkvinner.*skvinner;

levealder_menn = sum(dot(alder,pr_menn));
levealder_kvinner = sum(dot(alder,pr_kvinner));

figure;
semilogy(alder, qmenn, '-');
xlabel('Alder [år]')
ylabel('log_{10}(Dødssannsynlighet)')
hold on
semilogy(alder, qkvinner, '--');
legend('Menn','Kvinner')
hold off
figure;
plot(alder, qmenn, '-');
xlabel('Alder [år]')
ylabel('Dødssannsynlighet')
hold on
plot(alder, qkvinner, '--');
legend('Menn','Kvinner')
hold off

figure
plot(alder, smenn, '-');
xlabel('Alder [år]')
ylabel('Overlevelsessannsynlighet')
hold on
plot(alder, skvinner, '--');
legend('Menn','Kvinner')
hold off

figure
semilogy(alder, Pr_menn, '-');
xlabel('Alder [år]')
ylabel('log_{10}(Punktsannsynlighet)')
hold on
semilogy(alder, Pr_kvinner, '--');
legend('Menn','Kvinner')
hold off
figure
plot(alder, Pr_menn, '-');
xlabel('Alder [år]')
ylabel('Punktsannsynlighet')
hold on
plot(alder, Pr_kvinner, '--');
legend('Menn','Kvinner')
hold off