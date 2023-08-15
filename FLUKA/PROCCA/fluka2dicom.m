function [dose] = fluka2dicom(sourcepath, nbins)

%% Read .dat file
fileID = fopen(sourcepath, 'r');

formatSpec = '%f %f %f %f %f %f %f %f %f %f';
sizeA = [10 Inf];

A = fscanf(fileID, formatSpec, sizeA);
A = A.';

fclose(fileID);

%% Reshape Fortran matrix A into a 3-D dose matrix  
B = A.';
dose = reshape(B(1:end-8), [nbins.x nbins.y nbins.z]);
% dose = dose .* 1.602176462*10^(-7) * 10^9; % in nGy

end