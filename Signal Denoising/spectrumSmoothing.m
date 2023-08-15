function [SI, lambda] = spectrumSmoothing(path, minval, maxval)

%% 
[~, name, ~] = fileparts(path);

data = readXLSXdocument(path);
lambda = cell2mat(data(:,1));
SI = cell2mat(data(:,2));
Fs = 1./mean(abs(diff(lambda)));

% save('SI_and_lambda.mat', 'SI', 'lambda')

figure();
plot(lambda, SI, '-')
xlabel('Wavelength')
ylabel('Signal Intensity')
ylim([minval maxval])
title([name ' spectrum'])

%% Check for values less and larger than minvalue and maxvalue, respectively
% if SI(1) <= minval
%     SI(1) = minval;
% end

ind_min = find(SI <= minval);
ind_max = find(SI >= maxval);
lambda([ind_min.', ind_max.']) = [];
SI([ind_min.', ind_max.']) = [];
% SI = fillmissing(SI, 'pchip');

figure();
plot(lambda, SI, '-')
xlabel('Wavelength')
ylabel('Signal Intensity')
ylim([minval maxval])
title([name ' spectrum'])

%% Define threshold difference between signal values
diffSI = abs(diff(SI));
thres = quantile(diffSI, 0.25);

%% Remove all sharp peaks

% % Remove all signal values differing from neighbouring values above 
% % computed threshold and replace their values estimated by 
% % shape-preserving piecewise cubic spline interpolation 
for k = 1:2
        for i = (1 + k):(length(SI) - k)
            if abs(SI(i) - SI(i+k)) > thres
                SI(i) = NaN;
                SI(i+k) = NaN;
            elseif abs(SI(i) - SI(i-k)) > thres
                SI(i) = NaN;
                SI(i-k) = NaN;
            end
        end
        SI = fillmissing(SI, 'pchip'); 
end
figure()
plot(lambda, SI, '-')
xlabel('Wavelength')
ylabel('Signal Intensity')
ylim([minval maxval])
title('Signal after interpolation')
title([name ' spectrum'])

%%
windowWidth = 25;
SI = movmedian(SI, windowWidth); % movemean movmin
%SI = imopen(SI, ones(25,1));
figure();
plot(lambda, SI, '-')
xlabel('Wavelength')
ylabel('Signal Intensity')
ylim([minval maxval])
title([name ' spectrum'])

%%
windowWidthSG = 47; % 55
SI = sgolayfilt(SI, 2, windowWidthSG);
figure();
plot(lambda, SI, '-')
xlabel('Wavelength')
ylabel('Signal Intensity')
ylim([minval maxval])
title([name ' spectrum'])

%% Envelope the signal by lower and upper peak traces and compute mean signal
% [envHigh, envLow] = envelope(SI, 22, 'peak');
% SI_mean = (envHigh + envLow)./2;
% figure()
% plot(lambda, SI, '-.', lambda, SI_mean, lambda, envHigh, lambda, envLow)
% xlabel('Wavelength')
% ylabel('Signal Intensity')
% ylim([minval maxval])
% title([name ' spectrum'])
% legend('SI', 'Mean', 'High', 'Low', 'location', 'best')
  
% figure()
% plot(lambda, SI_mean, '-')
% xlabel('Wavelength')
% ylabel('Signal Intensity')
% ylim([minval maxval])
% title([name ' spectrum'])

end