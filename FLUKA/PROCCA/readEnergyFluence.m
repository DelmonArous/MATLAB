function [E, dPhi] = readEnergyFluence(path, ndetectors, nbins)

%% Read relevant data; differential fluence vs. energy

data = readXLSXdocument(path);

E = [];
dPhi = [];
skip = 2;   % skip empty space

for i = 0:(ndetectors-1)
    
   E_start(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 1)) .* 10^6;     % in kV
   E_end(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 2)) .* 10^6;     % in kV
   dPhi(:, i+1) = cell2mat(data(3 + (nbins+skip+2)*i: ...
       3 + (nbins+skip+2)*i + nbins-1, 3));             % in 1/cm2/kV
   
   E(:, i+1) = (E_start(:,i+1) + E_end(:,i+1)) ./ 2; 
   
end

end