clear all;
close all;
clc;
% reads the audio file and store sound samples to the array S
[S fs] = wavread('castanets.wav');
x_s = S(:,1); % only the first audio channel
x_s = x_s./max(abs(x_s)); % ensures that the samples are between -1 and 1
x = x_s(1:2^17); % the 2 ^ 17 first samples
N = length(x);
y = fft(x); % FFT of x
threshold = 4000;
for n = 1:(N/2+1) % NB: Matlab indexing
    % relation between the DFT index and frequency
    % if ">", removes frequencies higher than the threshold
    % if "<", remove frequencies less than the threshold
    if (n*fs/N < threshold)
        y(n) = 0; % removes frequencies with index n
        y(N-n) = 0; % remove frequencies associated with N-n due to aliasing
    end
end
xnew = real(ifft(y)); % Inverse FFT of y
xnew = xnew/max(abs(xnew));
% Play the backtransformed signal, xnew
playerobj = audioplayer(xnew, fs);
playblocking(playerobj);
% Plot
plot(1:N, abs(y))
xlabel('n')
ylabel('|y_n|')