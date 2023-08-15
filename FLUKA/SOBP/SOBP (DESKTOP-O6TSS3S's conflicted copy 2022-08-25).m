clear all;
close all;
fclose('all');
clc;

%% Path
destpath = 'C:\Users\delmo\Desktop\temp';

%% Parameters
E0      = 0.027;    % maximum energy in GeV
alpha   = 0.022;
p0      = 1.77;
n       = 20;       % number of intervals in SOBP
khi     = 0.75;     % fraction of the full range taken by the SOBP

%% Full range
R0 = R(E0, alpha, p0);

%% Initialization
r = [];
e = [];
w = [];

%% Loop through proton energies

for k = 0:n
    
    r(k+1) = (1 - (1 - (k/n)) * khi) * R0; 
    e(k+1) = (r(k+1)/alpha)^(1/p0);
    
    if k == 0
        w(k+1) = 1 - (1 - (1/(2*n)))^(1 - (1/p0));
    elseif k == n
        w(k+1) = (1/(2*n))^(1 - (1/p0));
    else
        w(k+1) = (1 - (1/n)*(k - 0.5))^(1 - (1/p0)) - ...
            (1 - (1/n)*(k + 0.5))^(1 - (1/p0));
    end
    
end

%% Define energy interval for each proton beamlet
e_min = [e(1)-mean(diff(e)) e(1:end-1)];
e_max = e;

%% Write to file
fid = fopen(fullfile(destpath, ['SOBP_khi' num2str(khi*100) '.txt']), 'wt');
fprintf(fid, '%.6f  %.6f   %.6f\n', [e_min; e_max; w]);
fclose(fid);
