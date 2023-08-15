function [E, dPhi, dPhi_fit, mu] = readEnergyFluence(sourcepath, ndetectors, nbins)

%% Read relevant data; differential fluence vs. energy

data = readXLSXdocument(sourcepath)

E = [];
dPhi = [];
skip = 2;   % skip empty space

for i = 0:(ndetectors-1)
    
   E_start(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 1)) .* 1000;     % in MeV
   E_end(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 2)) .* 1000;     % in MeV
   dPhi(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 3));             % in 1/cm2/GeV
   
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

end