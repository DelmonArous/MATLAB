function [mu] = plotEnergyFluence(sourcepath, filename, str, ndetectors, nbins)

%% Read relevant data; differential fluence vs. energy

data = readXLSXdocument(fullfile(sourcepath, filename));

E = [];
dPhi = [];
skip = 2;   % skip empty space

for i = 0:(ndetectors-1)
    
   E_start(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 1)) .* 10^3;     % in MeV
   E_end(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 2)) .* 10^3;     % in MeV
   dPhi(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 3));             % in 1/cm2/MeV
   
   E(:, i+1) = (E_start(:,i+1) + E_end(:,i+1)) ./ 2; 
   
end

%% Fit Gaussian distribution to data

f_gauss1 = {};
dPhi_fit = [];
mu = [];

for i = 1:ndetectors
    
    f_gauss1{i} = fit(E(:,i), dPhi(:,i), 'gauss1');
    dPhi_fit(:,i) = f_gauss1{i}.a1 .* exp(-((E(:,i) - ...
        f_gauss1{i}.b1)./f_gauss1{i}.c1).^2);
    mu = [mu f_gauss1{i}.b1];
    
end

%% Plot
lineSpec = {'-r', '-g', '-b', '-c', '-m', '-k', '-y'};
figure();
hold on
for i = 1:ndetectors
    plot(E(:,i), dPhi(:,i), lineSpec{i}, 'LineWidth', 1.0)    
    str{i} = [str{i} ' (' num2str(round(mu(i),2)) ' MeV)'];
end
hold off
ylim([0 inf])
xlabel('Energy (MeV)')
ylabel('Relative fluence (#/cm^2/MeV)') % Relative fluence (1/cm^2/kV)
lgd = legend(str, 'Location', 'NorthWest');
title(lgd, '15.2 MeV proton')
set(gca, 'FontSize', 14)


end