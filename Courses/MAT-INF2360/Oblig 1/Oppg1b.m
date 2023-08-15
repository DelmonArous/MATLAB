clear all;
clc;
% reads the audio file and store sound samples to the array S
[S fs] = wavread('castanets.wav');
x_s = S(:,1); % only the first audio channel, ie first column
x_s = x_s./max(abs(x_s)); % ensures that the samples are between -1 and 1
% initializing arrays
stop_DFTImpl = [];
stop_FFTImpl = [];
stop_fft = [];
n_ = [];
% loop over the 2^n first samples
for n = 4:14
    x = x_s(1:2^n); % get the 2^n samples
    n_ = [n_ n];
    % DFT implementation
    if (n <= 12) % to avoid error due to full system memory
        start_DFTImpl = tic;
        y_DFTImpl = DFTImpl(x);
        stop_DFTImpl = [stop_DFTImpl toc(start_DFTImpl)];
    end
    % FFT implementation
    start_FFTImpl = tic;
    y_FFTImpl = FFTImpl(x);
    stop_FFTImpl = [stop_FFTImpl toc(start_FFTImpl)];
    % Matlab's built in fft
    start_fft = tic;
    y_fft = fft(x);
    stop_fft = [stop_fft toc(start_fft)];
end
% Plot
plot(n_(1:9), stop_DFTImpl, 'o-' , n_, stop_FFTImpl, 'x-', n_, stop_fft, '*-')
xlabel('n')
ylabel('Tidsforbruk [s]')
legend('DFTImpl(x)', 'FFTImpl(x)', 'fft(x)')