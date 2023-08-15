clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source directory of files

path_MC_dose = { ... 
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v6\X-ray 220 kV spectrum wGRID\Dose maps\DAT - total dose', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v6\X-ray 220 kV spectrum wGRID\Dose maps\DAT - e dose', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v6\X-ray 220 kV spectrum wGRID\Dose maps\DAT - EM dose'};
path_MC_fluence = { ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v6\X-ray 220 kV spectrum wGRID\Fluence maps\DAT - gamma fluence', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\ProGRID\X-ray Results\v6\X-ray 220 kV spectrum wGRID\Fluence maps\DAT - e fluence'};

%% Variables

pattern = {'open', 'GRID'};

%% Read dose maps

for i = 2 % 1:length(path_MC_dose)
    
    filelist        = getAllFiles(path_MC_dose{i});
    dosemap_sum     = zeros(733,508); 
    flask_count     = 0;
    
    for j = 1:length(filelist)
        
        [filepath, filename, ext] = fileparts(filelist{j});
        [DoseMap_MC, RI_MC] = read2DMap(filelist{j});
        DoseMap_MC = DoseMap_MC .* 1.602176462*10^(-7); % * 10^9; % in nGy
        
        % Flip dose maps relative to flask 1
        if contains(filename, "Flask1", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Dose = flip(flip(DoseMap_MC, 1), 2);
        elseif contains(filename, "Flask2", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Dose = flip(DoseMap_MC, 1);
        elseif contains(filename, "Flask4", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Dose = flip(DoseMap_MC, 2);
        else
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Dose = DoseMap_MC;
        end
        
        MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', j)). ...
            RI = RI_MC;
        MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', j)). ...
            filename = filename;
        
        dosemap_sum = dosemap_sum + ...
            MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', j)).Dose;
        flask_count = flask_count + 1;
        
    end
    
    % Average dose map over all four flasks
    MC.(sprintf('%s', pattern{2})).MeanDoseMap = dosemap_sum./flask_count;
    
    % Spatial resolution in cm/pixel
    RI_MC       = MC.(sprintf('%s', pattern{2})).Flask3.RI;
    pxsize_MC   = [ diff(RI_MC.YWorldLimits) / size(MC.(sprintf('%s', pattern{2})).MeanDoseMap,1) ...
        diff(RI_MC.XWorldLimits) / size(MC.(sprintf('%s', pattern{2})).MeanDoseMap,2) ];
    
    % Plot
    h = figure();
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(MC.(sprintf('%s', pattern{2})).MeanDoseMap, ...
        MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', 3)).RI, ...
        [], 'colormap', jet(256))
    c = colorbar;
    c.Label.String = 'Dose (Gy per unit primary weight)';
    xlabel('cm'); ylabel('cm');
    grid on
    set(gca, 'FontSize', 16)
    shading interp
    
end

%% Read fluence maps

for i = 2 % 1:length(path_MC_fluence)
    
    filelist        = getAllFiles(path_MC_fluence{i});
    fluencemap_sum  = zeros(733,508); 
    flask_count     = 0;
    
    for j = 1:length(filelist)
        
        [filepath, filename, ext] = fileparts(filelist{j});
        [FluenceMap_MC, RI_MC] = read2DMap(filelist{j});
        
        % Flip fluence maps relative to flask 1
        if contains(filename, "Flask1", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Fluence = flip(flip(FluenceMap_MC, 1), 2);
        elseif contains(filename, "Flask2", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Fluence = flip(FluenceMap_MC, 1);
        elseif contains(filename, "Flask4", 'IgnoreCase', true)
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Fluence = flip(FluenceMap_MC, 2);
        else
            MC.(sprintf('%s', pattern{2})). ...
                (sprintf('Flask%i', j)).Fluence = FluenceMap_MC;
        end
        
        MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', j)). ...
            RI = RI_MC;
        MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', j)). ...
            filename = filename;
        
        fluencemap_sum = fluencemap_sum + ...
            MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', j)).Fluence;
        flask_count = flask_count + 1;
        
    end
    
    % Average dose map over all four flasks
    MC.(sprintf('%s', pattern{2})).MeanFluenceMap = fluencemap_sum./flask_count;
    
    % Spatial resolution in cm/pixel
    RI_MC       = MC.(sprintf('%s', pattern{2})).Flask3.RI;
    pxsize_MC   = [ diff(RI_MC.YWorldLimits) / size(MC.(sprintf('%s', pattern{2})).MeanFluenceMap,1) ...
        diff(RI_MC.XWorldLimits) / size(MC.(sprintf('%s', pattern{2})).MeanFluenceMap,2) ];
    
    % Plot
    h = figure();
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(MC.(sprintf('%s', pattern{2})).MeanFluenceMap, ...
        MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', 3)).RI, ...
        [], 'colormap', jet(256))
    c = colorbar;
    c.Label.String = 'Fluence (particles/cm^2 per unit primary weight)';
    xlabel('cm'); ylabel('cm');
    grid on
    set(gca, 'FontSize', 16)
    shading interp
    
end

%% LET

LETmap = 1.33*((6.242*10^(15))/10000000) .* ...
    (MC.GRID.MeanDoseMap./MC.GRID.MeanFluenceMap);
LETmap(isnan(LETmap)) = 0;

% Plot
h = figure();
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
imshow(LETmap, ...
    MC.(sprintf('%s', pattern{2})).(sprintf('Flask%i', 3)).RI, ...
    [0 40], 'colormap', jet(256))
c = colorbar;
c.Label.String = 'LET (keV/\mum)';
xlabel('cm'); ylabel('cm');
grid on
set(gca, 'FontSize', 16)
shading interp
