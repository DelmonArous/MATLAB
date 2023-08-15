function [b, bint, stats] = EstimateConcentrationsBeerLambertLaw(sourcefile)

%% Read data
[~, n, e] = fileparts(sourcefile);
data = readXLSXdocument(sourcefile);

% Check if there are at lest 2 absorbtion spectra
if size(data,2) < 9
   error(['Error: The data file ' n e ' does not contain enough spectra']);
end

n_spectra = (size(data,2) - 3) / 3; % nr. of spectra
n_values = size(data,1);
spectra_names = {};
lambda_array = [];
D_array = [];
lineSpec = {'r', 'g', 'b', 'c', 'm', 'k', 'y'};

for i = 1:(n_spectra+1)
    
    spectra_names{i} = cell2mat(data(1, 3*i-2));
    lambda_array(:,i) = cell2mat(data(1:end, 3*i-1));
    D_array(:,i) = cell2mat(data(1:end, 3*i));
    
end

%% Check if data contains NaNs 
temp_lambda = [];
temp_D = [];
NaNcolcounter = [];
nonNaNcolcounter = [];

if sum(isnan(lambda_array(:))) > 0
    NaNcolcounter = 0;
    nonNaNcolcounter = 0;
    for i = 1:n_spectra
        
        if sum(isnan(lambda_array(:,i))) == 0
            NaNcolcounter = NaNcolcounter + 1;
            temp_lambda(:,NaNcolcounter) = lambda_array(:,i);
            temp_D(:,NaNcolcounter) = D_array(:,i);
            NaNcolumns(NaNcolcounter) = i;
        else
            nonNaNcolcounter = nonNaNcolcounter + 1;
            nonNaNcolumns(nonNaNcolcounter) = i;
        end
        
    end
    
    lambda_array = rmmissing(lambda_array);
    D_array = rmmissing(D_array);
    
end

%% 1D interpolation
for i = 1:length(NaNcolumns)
    
    temp_D = interp1(temp_lambda(:,i), temp_D(:,i), lambda_array(:,nonNaNcolumns(1)), ...
        'pchip', 'extrap');
    D_array(:, NaNcolumns(i)) = temp_D;
    lambda_array(:, NaNcolumns(i)) = lambda_array(:,nonNaNcolumns(1));
    
end

%% Check if the finalized wavelength columns are equal
for i = 1:size(lambda_array,1)
    
    if any(diff(lambda_array(i,:)))
            error(['The spectra wavelengths in the data file ' n e ' are unequal.'])
    end

end

%% Check if the finalized array is correct
% lambda_PpIX = rmmissing(cell2mat(data(1:end, 2)));
% D_PpIX = rmmissing(cell2mat(data(1:end, 3)));
% lambda_Ppp = rmmissing(cell2mat(data(1:end, 5)));
% D_Ppp = rmmissing(cell2mat(data(1:end, 6)));
% lambda_tot = rmmissing(cell2mat(data(1:end, 8)));
% D_tot = rmmissing(cell2mat(data(1:end, 9)));
% % 
% D_Ppp = interp1(lambda_Ppp, D_Ppp, lambda_PpIX, 'pchip', 'extrap');
% lambda_Ppp = lambda_PpIX;
% D_smh = [D_PpIX D_Ppp D_tot];
% lambda_smh = [lambda_PpIX lambda_Ppp lambda_tot];
% D_smh == D_array 
% lambda_smh == lambda_array

%% Design matrix
X = ones(size(D_array(:,1)));  % X = [ones(size(D_PpIX)) D_PpIX D_Ppp];
for i = 1:n_spectra
    X = [X D_array(:,i)];
end

%% Test and multiple regression (assume D_i = epsilon_i)
c = [1 1]; % 1.02 2.5 in 10^(-6) M
D_test = zeros(size(X,1), 1);
for i = 1:n_spectra
   D_test = D_test + c(i) .* X(:,i+1); 
end

[b, bint, ~, ~, stats] = regress(D_array(:,end), X); % Removes NaN data % D_array(:,end)
Y = X*b;
maxerr = max(abs(Y - D_array(:,end))); % max absolute value of the deviation % D_array(:,end)

%% Plot
figure()
hold on
for i = 1:(n_spectra)

    plot(lambda_array(:,i), D_array(:,i), lineSpec{i})

end
plot(lambda_array(:,1), D_test, lineSpec{3})
hold off
ylim([0 900])
xlabel('wavelength (nm)')
ylabel('Absorbance (a.u.)')
legend(spectra_names, 'Location', 'best');

figure()
hold on
for i = 1:(n_spectra+1)

    plot(lambda_array(:,i), D_array(:,i), lineSpec{i})

end
hold off
ylim([0 900])
xlabel('wavelength (nm)')
ylabel('Absorbance (a.u.)')
legend(spectra_names, 'Location', 'best');


end