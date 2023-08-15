x = [3.9 5.6 4.1 4.2 4.0 3.6 5.9 4.5 3.6 5.0 2.9 4.3];

gjennomsnitt = sum(x)/length(x);
empirisk_median = median(x);
empirisk_standardavvik = std(x);
kvartildifferanse = iqr(x);

edges = [2, 3, 4, 5, 6];
k = histc(x,edges);
k = normpdf(k);
bar(edges, k, 'histc');
figure(gcf);

figure;
hist(x,sqrt(length(x)));
xlabel('Vitalkapasitet [liter]');
ylabel('Antall personer med tilhørende vitalkapasitet');