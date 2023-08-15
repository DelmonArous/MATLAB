function [depth, dose] = getDepthDose(sourcefile)

% Read depth dose data
A = importdata(sourcefile);

% Initialize
bin_start   = [];
bin_end     = [];
dose        = [];

% Compute depth in cm and dose in nGy
bin_start   = A.data(1:end, 1);
bin_end     = A.data(1:end, 2);
dose        = A.data(1:end, 3) .* 1.602176462*10^(-7) .* 10^9;  % in nGy;
depth       = (bin_start + bin_end) ./ 2;                       % in cm

end