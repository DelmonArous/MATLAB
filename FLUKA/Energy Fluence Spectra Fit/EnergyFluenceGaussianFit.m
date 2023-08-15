function [] = EnergyFluenceGaussianFit(path, sourcefilename, destfilename, ...
    ndetectors, nbins, str)

%% Read relevant data; differential fluence vs. energy

data = readXLSXdocument(fullfile(path,sourcefilename));


ndetectors = 2;
nbins = 400;
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


%% Plot and fit Gaussian distribution to data

% gauss1Eqn = fittype('a*exp(-((x - b)/(2*c))^2)', ...
%     'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a','b','c'});
% gauss2Eqn = fittype('a1*exp(-((x - b1)/(2*c1))^2) + a2*exp(-((x - b2)/(2*c2))^2)', ...
%     'dependent', {'y'}, 'independent', {'x'}, ...
%     'coefficients', {'a1','b1','c1','a2','b2','c2'});
% gauss1Eqn = 'a*exp(-0.5*((x - b)/c)^2)';
% gauss2Eqn = 'a1*exp(-0.5*((x - b1)/c1)^2) + a2*exp(-0.5*((x - b2)/c2)^2)';
% startPoints1 = {[1.975 12.54 0.1802], [0.3692 3.586 0.5957], ...
%     [0.316 2.923 0.6921], [0.2508 2.143 0.8498], [0.192 1.495 0.9537]};
% startPoints2 = {[2.373 12.58 0.2023 -0.7384 12.7 0.1679], ...
%     [0.2559 3.709 0.4978 0.1512 3.347 0.6217], ...
%     [0.2307 3.065 0.5659 0.1228 2.604 0.738], ...
%     [0.1006 1.697 0.9046 0.1901 2.335 0.6598], ...
%     [0.1955 1.63 1.114 -0.04926 2.635 0.6956]};

f_gauss1 = {};
f_gauss2 = {};
temp_dPhi_gauss1 = [];
temp_dPhi_gauss2 = [];
mu = []; SD = [];

lineSpecData = {'ro', 'go', 'bo', 'co', 'mo', 'ko', 'yo'};
lineSpecFit = {'-r', '-g', '-b', '-c', '-m', '-k', '-y'};
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
figure();

for i = 1:ndetectors
    
    f_gauss1{i} = fit(E(:,i), dPhi(:,i), 'gauss1');
    f_gauss2{i} = fit(E(:,i), dPhi(:,i), 'gauss2');     
    
    temp_dPhi_gauss1(:,i) = f_gauss1{i}.a1 .* exp(-((E(:,i) - f_gauss1{i}.b1)./f_gauss1{i}.c1).^2);
    temp_dPhi_gauss2(:,i) = f_gauss2{i}.a1 .* exp(-((E(:,i) - f_gauss2{i}.b1)./f_gauss2{i}.c1).^2) + ...
        f_gauss2{i}.a2 .* exp(-((E(:,i) - f_gauss2{i}.b2)./f_gauss2{i}.c2).^2);
    
    hold on
    %plot(E(:,i), dPhi(:,i), lineSpecData{i}, 'MarkerSize', 1);
    plot(E(:,i), temp_dPhi_gauss1(:,i), lineSpecFit{i}); % p(2*i)
    %plot(f_gauss1{i}, lineSpecFit{i})

    mu = [mu; (f_gauss2{i}.b1 + f_gauss2{i}.b2)/2];
    SD = [SD; sqrt(((f_gauss2{i}.c1/2)^2 + (f_gauss2{i}.c2/2)^2)/2)];
    
    % Append text to legend string
    lgdstr{i} = [str{i} ' (' num2str(round(f_gauss1{i}.b1,1)) ' MeV)']; 
    
    % Append text to .xlsx header string
    header{2*i-1} = 'Energy (MeV)';
    xlsx_col{2*i-1} = [alphabet(2*i-1) '2'];
    header{2*i} = ['Relative fluence (1/cm2/GeV): ' str{i}];
    xlsx_col{2*i} = [alphabet(2*i) '2'];

end
legend(lgdstr, 'Location', 'best') % p(2:2:end)
xlabel('Energy (MeV)')
ylabel('Relative fluence (1/cm^2/GeV)')
ylim([0 inf])

%% Write to XLSX file

xlswrite(fullfile(path, destfilename), header, 1, 'A1')
for i = 1:ndetectors
    xlswrite(fullfile(path, destfilename), E(:,i), 1, xlsx_col{2*i-1})
    xlswrite(fullfile(path, destfilename), temp_dPhi_gauss2(:,i), 1, xlsx_col{2*i})
end

f_gauss2{:}

end