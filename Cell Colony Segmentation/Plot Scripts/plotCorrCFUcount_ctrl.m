clc
clear all
close all
warning('off','all')

sourcepath = 'C:\Users\delmo\Desktop\CountCTRL\T47D';
folderList = getAllFolders(sourcepath);

%% Store all count data in a structure
struct = {};
for i = 1:length(folderList)
    
    [~, foldername, ~] = fileparts(folderList{i});
    fileList = getAllFiles(folderList{i});
    
    for j = 1:length(fileList)
        
        [path, filename, ext] = fileparts(fileList{j});
        data = readXLSXdocument(fileList{j});
        struct(i).(sprintf('%s', strrep(string(filename), '-', ''))).Count ...
            = cell2mat(data(2:end, 2));
        
    end
    
end

if isempty(struct)
    prcnt_diff('Error occured. No colony data were found.');
end

% % % Dose points in Gy
% % Ingunn T47D: [0 0.1 0.3 0.5 0.75 1 2 5]
% % Hilde A549: [0 0.5 1 2 5 10]
% % doses = [0 0.5 1 2 5 10];
% % dose = repelem(doses, 4);
% % dose = [0 0 0 0 dose].';
% % Ingunn T47D: [repelem(200,length(dose)-4) repelem(400,4)]
% % Hilde A549: [repelem(200,length(dose)-8) repelem(800,4) repelem(20000,4)]
% % cells_seeded = [repelem(200,length(dose)-8) repelem(800,4) repelem(20000,4)];
% % n_ctrl = 8;
% % Ingunn T47D: [1.215686 1.236364]
% % Hilde A549: [301/229]
% % M = [301/229]; % multiplicity

n_ctrl = 8;

%% Correlate automatic and manual count
for i = 1:length(struct)
    
    CFUcount = [];
    
    for field = fieldnames(struct(i))'
        if ~isempty(struct(i).(field{1}))
            (field{1})
            % Automatic and manual count are stored in the 1st and 2nd
            % column, respectively
            CFUcount = [CFUcount struct(i).(field{1}).Count];
        end
    end
    
    auto = CFUcount(:,1);
    manual = CFUcount(:,2);
    
%     autoMean = [mean(auto(1:n_ctrl)) auto(n_ctrl+1:end).'];
%     autoMean = [autoMean(1) arrayfun(@(i)mean(auto(i:i+4-1)), ...
%         (n_ctrl+1):4:length(auto)-4+1)];
%     manualMean = [mean(manual(1:n_ctrl)) manual((n_ctrl+1):end).'];
%     manualMean = [manualMean(1) arrayfun(@(i)mean(manual(i:i+4-1)), ...
%         (n_ctrl+1):4:length(manual)-4+1)];
%     autoStd = [std(auto(1:n_ctrl)) auto((n_ctrl+1):end).'];
%     autoStd = [autoStd(1) arrayfun(@(i)std(auto(i:i+4-1)), ...
%         (n_ctrl+1):4:length(auto)-4+1)];
%     manualStd = [std(manual(1:n_ctrl)) manual((n_ctrl+1):end).'];
%     manualStd = [manualStd(1) arrayfun(@(i)std(manual(i:i+4-1)), ...
%         (n_ctrl+1):4:length(manual)-4+1)];
    
%     % Estimate plating efficiency (PE)
%     PE_auto = autoMean(1)/cells_seeded(1);
%     PE_manual = manualMean(1)/cells_seeded(1);
%     
%     % Estimate individual surviving fractions (SFs)
%     SF_auto = auto((n_ctrl+1):end) ./ ...
%         (cells_seeded((n_ctrl+1):end)' .* PE_auto);
%     SF_auto = [repelem(1,n_ctrl) SF_auto']';
%     SF_manual = manual((n_ctrl+1):end) ./ ...
%         (cells_seeded((n_ctrl+1):end)' .* PE_manual);
%     SF_manual = [repelem(1,n_ctrl) SF_manual']';
% 
%     SF_auto = (M(i) - sqrt(M(i)^2 - 4.*(M(i)-1).*SF_auto)) ./ (2*(M(i)-1));
%     SF_manual = (M(i) - sqrt(M(i)^2 - 4.*(M(i)-1).*SF_manual)) ./ (2*(M(i)-1));
    
    %% Linear correlation
    [r, pval] = corr(auto, manual, 'type', 'Pearson');
    str = strcat('n=', num2str(length(auto)), ...
        ', R^2=', num2str(round(r^2,3)), ...
        ', P=', num2str(round(pval,9)));
    X = [ones(length(auto), 1) manual];
    beta = X\auto;
    pval
    
    % Plot
    figure();
    hold on
    plot(manual, auto, 'ko')
    plot(manual, X*beta, 'k-');
    hold off
    xlabel('Manual Count')
    ylabel('Automatic Count')
    lgd = legend(['T-47D data'], 'Location', 'NorthWest');
    %%%%% num2str(10^(floor(log(abs(pval_Q2))./log(10)) + 1)));
    title(lgd, str)
    set(gca, 'FontSize', 14)
    
%     figure();
%     hold on
%     errorbar(doses, autoMean, autoStd, '-or')
%     errorbar(doses, manualMean, manualStd, '-xb')
%     hold off
%     xlabel('Dose (Gy)')
%     ylabel('Colony Count')
%     lgd = legend('Automatic', 'Manual');
%     title(lgd, ['A549, Exp ' num2str(i)])
%     xlim([-0.1 max(doses)+0.1])
%     ylim([-0.5 Inf])
%     set(gca, 'FontSize', 14)
     
%     prcnt_diff = ((auto-manual)./((auto+manual)./2)) .* 100;
%     [curve, ~, ~] = fit(dose, prcnt_diff, 'smoothingspline');
%     figure();
%     hold on
%     h1 = plot(curve, dose, prcnt_diff);
%     plot(dose, repelem(10, length(dose)), '--k')
%     plot(dose, repelem(0, length(dose)), '-k')
%     plot(dose, repelem(-10, length(dose)), '--k')
%     %     scatter(dose, manual-auto, 'o', 'filled', 'MarkerEdgeColor', 'black')
%     %     plot(dose, f, '-k')
%     hold off
%     set(h1, 'MarkerSize', 16, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
%     lgd = legend('Data', 'Smoothing spline', 'Location', 'best');
%     title(lgd, ['T-47D, Exp ' num2str(i)])
%     xlabel('Dose (Gy)')
%     ylabel('% Difference (auto - manual)')
%     xlim([-0.1 max(doses)+0.1])
% %     ylim([min(prcnt_diff)-1 max(prcnt_diff)+1])
%     set(gca, 'FontSize', 14)    
   
      
end