function [pat] = getFiles(patients, modality, seriesdesc1, seriesdesc2, seriesdesc3)

% The function creates a structure of relevant patients with DICOM images
% taken with specified modality and series

% The following variables are required for proper execution:
%   patients: loadPatients object of allocated patients
%   modality: string containing the modality of the DICOM image to be
%               searched for and allocated
%   series: string containing the series of the DICOM image to be searched
%               for and allocated

counter = 0;
for i = 1:length(patients)
    
    pat{i}.name = patients{i}.name;
    
    for j = 1:length(patients{i}.img)
        
        if nargin == 2
            
            if (isfield(patients{i}.img{j}, 'Modality') && ...
                    strcmpi(patients{i}.img{j}.Modality, modality))
                counter = counter + 1;
                pat{i}.img{counter} = patients{i}.img{j};
            end
            
        elseif nargin == 3
            
            if (isfield(patients{i}.img{j}, 'Modality') && ...
                    strcmpi(patients{i}.img{j}.Modality, modality) && ...
                    isfield(patients{i}.img{j}, 'SeriesDescription') && ...
                    contains(patients{i}.img{j}.SeriesDescription, seriesdesc1, 'IgnoreCase', true))
                counter = counter + 1;
                pat{i}.img{counter} = patients{i}.img{j};
            end
            
        elseif nargin == 4
            
            if ( isfield(patients{i}.img{j}, 'Modality') && ...
                    strcmpi(patients{i}.img{j}.Modality, modality) && ...
                    isfield(patients{i}.img{j}, 'SeriesDescription') && ...
                    contains(patients{i}.img{j}.SeriesDescription, seriesdesc1, 'IgnoreCase', true) && ...
                    contains(patients{i}.img{j}.SeriesDescription, seriesdesc2, 'IgnoreCase', true) )
                counter = counter + 1;
                pat{i}.img{counter} = patients{i}.img{j};
            end
            
        elseif nargin == 5
            
            if ( isfield(patients{i}.img{j}, 'Modality') && ...
                    strcmpi(patients{i}.img{j}.Modality, modality) && ...
                    isfield(patients{i}.img{j}, 'SeriesDescription') && ...
                    contains(patients{i}.img{j}.SeriesDescription, seriesdesc1, 'IgnoreCase', true) && ...
                    (contains(patients{i}.img{j}.SeriesDescription, seriesdesc2, 'IgnoreCase', true) || ...
                    contains(patients{i}.img{j}.SeriesDescription, seriesdesc3, 'IgnoreCase', true)) )
                counter = counter + 1;
                pat{i}.img{counter} = patients{i}.img{j};
            end
            
        else
            error('Incorrect input arguments!')
        end
        
    end
    
    counter = 0;
    
end

% Remove empty cells in the structure
pat = pat(~cellfun('isempty', pat));

clear counter

end