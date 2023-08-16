clear all;

filnavn= 'soldata.txt';
fileID = fopen(filnavn, 'r');
A = fscanf(fileID,'%f %f', [2,inf]);

[M,N] = size(A);
z = (2/N)*fft(A(2,:));
T = A(1,N)-A(1,1);
deltaF = 1.0/T; 

for i = 1:N
    f(i) = deltaF*(i-1);
end;

figure;
plot(A(1,:),A(2,:),'-b');
xlabel('År')
ylabel('Antall solflekker')
figure;
plot(f,abs(z),'.-k');
xlabel('Frekvens (år^-^1)');
ylabel('Antall solflekker');
axis([min(f) max(f)/2 min(A(2,:)) max(A(2,:))])