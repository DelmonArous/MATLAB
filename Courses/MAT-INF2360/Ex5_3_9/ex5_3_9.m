clear all;
close all;
clc;

[S fs] = wavread('castanets.wav');
x = S(:,1);
x = x./max(abs(x));
x = x(1:2^17);
N = length(x);

xnew = DWTHaarImpl(x,2);

figure;
plot((0:N-1), x)
figure;
plot((0:N-1), xnew)
