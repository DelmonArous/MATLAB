clc
clear all
close all
warning('off','all')

sourcepath = 'C:\Users\delmo\Desktop\Excel\A549\X-ray\17122020';
folderList = getAllFolders(sourcepath);

struct = {};

for i = 1:length(folderList)
    
    [~, foldername, ~] = fileparts(folderList{i});
    fileList = getAllFiles(folderList{i});
    
    for j = 1:length(fileList)
        
        [path, filename, ext] = fileparts(fileList{j});
        data = readXLSXdocument(fileList{j});
        
        struct.Count.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            length(cell2mat(data(2:end, 1)));
        struct.Area.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 1));
        struct.MedianRed.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 4));
        struct.MedianGreen.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 5));
        struct.MedianBlue.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 6));
        struct.MedianGray.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 7));
        struct.MedianPCA1.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 8));
        struct.MedianPCA2.(sprintf('%s', strrep(foldername, ' ', ''))). ...
            (sprintf('%s', strrep(filename, '-', '_'))) = ...
            cell2mat(data(2:end, 9));
        %         struct(j).(sprintf('%s', strrep(foldername, ' ', ''))).StdRed ...
        %             = cell2mat(data(2:end, 10));
        %         struct(j).(sprintf('%s', strrep(foldername, ' ', ''))).StdGreen ...
        %             = cell2mat(data(2:end, 11));
        %         struct(j).(sprintf('%s', strrep(foldername, ' ', ''))).StdBlue ...
        %             = cell2mat(data(2:end, 12));
        %         struct(j).(sprintf('%s', strrep(foldername, ' ', ''))).StdGray ...
        %             = cell2mat(data(2:end, 13));
        %         struct(j).(sprintf('%s', strrep(foldername, ' ', ''))).StdPCA1 ...
        %             = cell2mat(data(2:end, 14));
        %         struct(j).(sprintf('%s', strrep(foldername, ' ', ''))).StdPCA2 ...
        %             = cell2mat(data(2:end, 15));
        
    end
    
end

if isempty(struct)
    error('Error occured. No colony data were found.');
end

%% Rearrange 
structvalues_GRID = {};
for field1 = fieldnames(struct)'
            
        counter = 0;
        
        for field2 = fieldnames(struct.(field1{1}))'
            if ~strcmp(field2{1}, 'Open')
                for filename = fieldnames(struct.(field1{1}).(field2{1}))'
                    counter = counter + 1;
                    structvalues_GRID(counter).(field1{1}) = ...
                        struct.(field1{1}).(field2{1}).(filename{1});
                end  
            end
        end
            
end
structvalues_Open = {};
for field1 = fieldnames(struct)'
            
        counter = 0;
        
        for field2 = fieldnames(struct.(field1{1}))'
            if ~strcmp(field2{1}, 'GRIDStripes')
                for filename = fieldnames(struct.(field1{1}).(field2{1}))'
                    counter = counter + 1;
                    structvalues_Open(counter).(field1{1}) = ...
                        struct.(field1{1}).(field2{1}).(filename{1});
                end
            end
        end
            
end

%% Plot
ylabel_str = {'Colony Count', 'Colony Area (mm^2)', ...
    'CFU Red Intensity (a.u)', ...
    'CFU Green Intensity (a.u)', ...
    'CFU Blue Intensity (a.u)', ...
    'CFU Grayscale Intensity (a.u)', ...
    'CFU PCA1 Intensity (a.u)', ...
    'CFU PCA2 Intensity (a.u)'};

dose = [0 2 5 10]; % in Gy
n_ctrl = 4;
counter = 0;
for variable = fieldnames(structvalues_Open)' 

    counter = counter + 1;
    val_GRID = [];
    val_Open = [];
    
    for i = 1:length(structvalues_Open)
        val_Open = [val_Open mean(structvalues_Open(i).(variable{1}))];
        val_GRID = [val_GRID mean(structvalues_GRID(i).(variable{1}))];
    end    
    
    val_Open_mean = [mean(val_Open(1:n_ctrl)) val_Open(n_ctrl+1:end)];
    val_Open_mean = [val_Open_mean(1) arrayfun(@(i)mean(val_Open(i:i+4-1)), ...
        (n_ctrl+1):4:length(val_Open)-4+1)];    
    val_GRID_mean = [mean(val_GRID(1:n_ctrl)) val_GRID(n_ctrl+1:end)];
    val_GRID_mean = [val_GRID_mean(1) arrayfun(@(i)mean(val_GRID(i:i+4-1)), ...
        (n_ctrl+1):4:length(val_GRID)-4+1)];
    
    val_Open_std = [std(val_Open(1:n_ctrl)) val_Open((n_ctrl+1):end)];
    val_Open_std = [val_Open_std(1) arrayfun(@(i)std(val_Open(i:i+4-1)), ...
        (n_ctrl+1):4:length(val_Open)-4+1)];
    val_GRID_std = [std(val_GRID(1:n_ctrl)) val_GRID((n_ctrl+1):end)];
    val_GRID_std = [val_GRID_std(1) arrayfun(@(i)std(val_GRID(i:i+4-1)), ...
        (n_ctrl+1):4:length(val_GRID)-4+1)];
    
    % Normalize to control
%     val_Open_mean = val_Open_mean ./ val_Open_mean(1); % uncomment for count plot
%     val_GRID_mean = val_GRID_mean ./ val_GRID_mean(1); %uncomment for count plot

%     val_Open_std = val_Open_std ./ val_Open_mean(1); % uncomment for count plot
%     val_GRID_std = val_GRID_std ./ val_GRID_mean(1); %uncomment for count plot
    
    figure();
    hold on
    errorbar(dose, val_Open_mean, val_Open_std, '-or')
    errorbar(dose, val_GRID_mean, val_GRID_std, '-xb')
    hold off
    xlabel('Dose (Gy)')
    ylabel(ylabel_str{counter})
    xlim([-0.1 max(dose)+0.1])
%     ylim([min(val_mean(:))-0.1 max(val_mean(:))+0.1])
    lgd = legend('Open', 'GRID stripes', 'Location', 'SouthWest');
    title(lgd, ['A549, Exp ' num2str(n_ctrl-1)])
    set(gca, 'FontSize', 14)
    
end

% for field = fieldnames(struct(1).GRIDStripes)'
%     
%     if ~strcmp(field{1}, 'Filename')
%         
%         counter = counter + 1;
%         value_GRID = []; value_Open = [];
%         
%         % First, loop over the controls
%         for j = find(~cellfun(@isempty, {struct.Control}))
%             value_GRID = [value_GRID mean(struct(j).Control.(field{1}))];
%             value_Open = [value_Open mean(struct(j).Control.(field{1}))];
%         end
%         % Then loop over the experiments
%         for j = 1:numel(struct)
%             value_GRID = [value_GRID mean(struct(j).GRIDStripes.(field{1}))];
%             value_Open = [value_Open mean(struct(j).Open.(field{1}))];
%         end
%         
%         valueMean_GRID = arrayfun(@(i) ...
%             mean(value_GRID(i:i+4-1)),1:4:length(value_GRID)-4+1)'; % n = 4
%         valueStd_GRID = arrayfun(@(i) ...
%             std(value_GRID(i:i+4-1)),1:4:length(value_GRID)-4+1)';
%         valueMean_Open = arrayfun(@(i) ...
%             mean(value_Open(i:i+4-1)),1:4:length(value_Open)-4+1)';
%         valueStd_Open = arrayfun(@(i) ...
%             std(value_Open(i:i+4-1)),1:4:length(value_Open)-4+1)';
%         
%         % Normalize to controls
%         if ~strcmp(field{1}, 'Count')
%             valueMean_GRID = valueMean_GRID ./ valueMean_GRID(1);
%             valueMean_Open = valueMean_Open ./ valueMean_Open(1);
%         end
%         
%         figure();
%         hold on
%         errorbar(dose, valueMean_GRID, valueStd_GRID, '-or')
%         errorbar(dose, valueMean_Open, valueStd_Open, '-xb')
%         xlabel('Dose (Gy)')
%         ylabel(ylabel_str{counter})
%         legend('GRID Stripes', 'Open')
%         xlim([0 max(dose)+1])
%         set(gca, 'FontSize', 12)
%         
%     end
%     
% end

%% Plot KDE
% xlabel_str = {'Colony Area (mm^2)', 'Colony Red Intensity (a.u)', ...
%     'Colony Green Intensity (a.u)', 'Colony Blue Intensity (a.u)', ...
%     'Colony Grayscale Intensity (a.u)', 'Colony PCA1 Intensity (a.u)', ...
%     'Colony PCA2 Intensity (a.u)'};
% lgd_strGRID = {'Control', 'Control', '2 Gy, GRID stripes', ...
%     '2 Gy, GRID stripes', '5 Gy, GRID stripes', '5 Gy, GRID stripes', ...
%     '10 Gy, GRID stripes', '10 Gy, GRID stripes'};
% lgd_strOpen = {'Control', 'Control', '2 Gy, Open', '2 Gy, Open', ...
%     '5 Gy, Open', '5 Gy, Open', '10 Gy, Open', '10 Gy, Open'};
% linespec_GRID = {'-k', '--k', '-r', '--r', '-r', '--r', '-r', '--r'};
% linespec_Open = {'-k', '--k', '-b', '--b', '-b', '--b', '-b', '--b'};
% counter = 0;
% for field = fieldnames(struct(1).GRIDStripes)'
%
%     if ~strcmp(field{1}, 'Filename') && ~strcmp(field{1}, 'Count')
%
%         counter = counter + 1;
%         KDE_GRID = []; KDE_Open = [];
%         valuesmesh_GRID = []; valuesmesh_Open = [];
%
%         % First, loop over the controls
%         for j = find(~cellfun(@isempty, {struct.Control}))
%             [~, KDE, valuesmesh, ~] = kde(struct(j).Control.(field{1}), ...
%                 2^12, min(struct(j).Control.(field{1})), ...
%                 max(struct(j).Control.(field{1})));
%             KDE_GRID = [KDE_GRID KDE];
%             KDE_Open = [KDE_Open KDE];
%             valuesmesh_GRID = [valuesmesh_GRID valuesmesh.'];
%             valuesmesh_Open = [valuesmesh_Open valuesmesh.'];
%         end
%         % Then loop over the experiments
%         for j = 1:numel(struct)
%             [~, KDE, valuesmesh, ~] = kde(struct(j).GRIDStripes.(field{1}), ...
%                 2^12, min(struct(j).GRIDStripes.(field{1})), ...
%                 max(struct(j).GRIDStripes.(field{1})));
%             KDE_GRID = [KDE_GRID KDE];
%             valuesmesh_GRID = [valuesmesh_GRID valuesmesh.'];
%
%             [~, KDE, valuesmesh, ~] = kde(struct(j).Open.(field{1}), ...
%                 2^12, min(struct(j).Open.(field{1})), ...
%                 max(struct(j).Open.(field{1})));
%             KDE_Open = [KDE_Open KDE];
%             valuesmesh_Open = [valuesmesh_Open valuesmesh.'];
%         end
%
%         % Estimate KDE
%         for j = 3:2:size(KDE_GRID,2)
%             figure();
%             hold on
%             plot(valuesmesh_GRID(:,1), KDE_GRID(:,1), linespec_GRID{1}, ...
%                 valuesmesh_GRID(:,2), KDE_GRID(:,2), linespec_GRID{2})
%             plot(valuesmesh_GRID(:,j), KDE_GRID(:,j), linespec_GRID{j}, ...
%                 valuesmesh_GRID(:,j+1), KDE_GRID(:,j+1), linespec_GRID{j+1}, ...
%                 valuesmesh_Open(:,j), KDE_Open(:,j), linespec_Open{j}, ...
%                 valuesmesh_Open(:,j+1), KDE_Open(:,j+1), linespec_Open{j+1})
%             ylabel('Kernel Density Estimation')
%             xlabel(xlabel_str{counter})
%             legend('Control', 'Control', ...
%                 lgd_strGRID{j}, lgd_strGRID{j+1}, ...
%                 lgd_strOpen{j}, lgd_strOpen{j+1}, 'Location', 'best')
%             hold off
%         end
%
%     end
%
% end

%% Plot histogram
% xlabel_str = {'Colony Area (mm^2)', 'Colony Red Intensity (a.u)', ...
%     'Colony Green Intensity (a.u)', 'Colony Blue Intensity (a.u)', ...
%     'Colony Grayscale Intensity (a.u)', 'Colony PCA1 Intensity (a.u)', ...
%     'Colony PCA2 Intensity (a.u)'};
% lgd_strGRID = {'Control', 'Control', '2 Gy, GRID stripes', ...
%     '2 Gy, GRID stripes', '5 Gy, GRID stripes', '5 Gy, GRID stripes', ...
%     '10 Gy, GRID stripes', '10 Gy, GRID stripes'};
% lgd_strOpen = {'Control', 'Control', '2 Gy, Open', '2 Gy, Open', ...
%     '5 Gy, Open', '5 Gy, Open', '10 Gy, Open', '10 Gy, Open'};
% linespec_GRID = {'-k', '--k', '-r', '--r', '-r', '--r', '-r', '--r'};
% linespec_Open = {'-k', '--k', '-b', '--b', '-b', '--b', '-b', '--b'};
% counter = 0;
% for field = fieldnames(struct(1).GRIDStripes)'
%
%     if ~strcmp(field{1}, 'Filename') && ~strcmp(field{1}, 'Count')
%
%         counter = counter + 1;
%         KDE_GRID = []; KDE_Open = [];
%         valuesmesh_GRID = []; valuesmesh_Open = [];
%
%         % First, loop over the controls
%         for j = find(~cellfun(@isempty, {struct.Control}))
%             histvalues = computeHistogram(struct(j).Control.(field{1}));
%             KDE_GRID = [KDE_GRID histvalues(:,2)];
%             KDE_Open = [KDE_Open histvalues(:,2)];
%             valuesmesh_GRID = [valuesmesh_GRID histvalues(:,1)];
%             valuesmesh_Open = [valuesmesh_Open histvalues(:,1)];
%         end
%         % Then loop over the experiments
%         for j = 1:numel(struct)
%             histvalues = computeHistogram(struct(j).GRIDStripes.(field{1}));
%             KDE_GRID = [KDE_GRID histvalues(:,2)];
%             valuesmesh_GRID = [valuesmesh_GRID histvalues(:,1)];
%
%             histvalues = computeHistogram(struct(j).Open.(field{1}));
%             KDE_Open = [KDE_Open histvalues(:,2)];
%             valuesmesh_Open = [valuesmesh_Open histvalues(:,1)];
%         end
%
%         % Estimate KDE
%         for j = 3:2:size(KDE_GRID,2)
%             figure();
%             hold on
%             plot(valuesmesh_GRID(:,1), KDE_GRID(:,1), linespec_GRID{1}, ...
%                 valuesmesh_GRID(:,2), KDE_GRID(:,2), linespec_GRID{2})
%             plot(valuesmesh_GRID(:,j), KDE_GRID(:,j), linespec_GRID{j}, ...
%                 valuesmesh_GRID(:,j+1), KDE_GRID(:,j+1), linespec_GRID{j+1}, ...
%                 valuesmesh_Open(:,j), KDE_Open(:,j), linespec_Open{j}, ...
%                 valuesmesh_Open(:,j+1), KDE_Open(:,j+1), linespec_Open{j+1})
%             ylabel('Percentage Colonies (%)')
%             xlabel(xlabel_str{counter})
%             legend('Control', 'Control', ...
%                 lgd_strGRID{j}, lgd_strGRID{j+1}, ...
%                 lgd_strOpen{j}, lgd_strOpen{j+1}, 'Location', 'best')
%             ylim([0 100])
%             set(gca, 'FontSize', 12)
%             hold off
%         end
%
%     end
%
% end
