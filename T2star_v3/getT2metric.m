function [R2starPatients] = getT2metric(T2wFatpath, T2wImgsourcepath, R2starPatients)

folderList = getAllFolders(T2wImgsourcepath);

for i = 1:length(folderList)
    
    [folderpath, foldername, ~] = fileparts(folderList{i});
    
    if contains(foldername, 'Patient')
%         [R2starave, T2met, ADC, fBV, HF_DWI] = ...
%             mapROI(T2wFatpath, fullfile(folderpath, foldername), R2starPatients);
        R2starPatients = ...
            mapROI(T2wFatpath, fullfile(folderpath, foldername), R2starPatients);
              
%         if ~isnan(R2starave) && ~isnan(T2met) && ~isnan(ADC) && ~isnan(fBV) && ~isnan(HF_DWI)
%             HF_DWIvec = [HF_DWIvec; HF_DWI]; 
%             ADCvec = [ADCvec; ADC];
%             fBVvec = [fBVvec; fBV];
%             R2starMedian = [R2starMedian; R2starave];
%             T2Metric = [T2Metric; T2met];
%         end
        
    end
    
end

% %% Plot R2* vs T2metric
% [r, pval] = corr(R2starMedian, T2Metric, 'type', 'Spearman');
% str = strcat('n=', num2str(length(T2Metric)), ...
%     ', r=', num2str(round(r,3)), ...
%     ', P=', num2str(round(pval,6)));
% 
% figure()
% hold on
% scatter(T2Metric, R2starMedian, [], fBVvec, 'filled');
% %plot(T2Metric, R2starMedian, 'bo')
% ylabel('R2* (s^{-1})')
% xlabel('T2 metric')
% lgd = legend('Tumor', 'Location', 'NorthEast');
% title(lgd, str)
% c = colorbar;
% c.Label.String = 'fBV';
% shading interp;
% jet(256);
% set(gca, 'FontSize', 11)
% 
% figure()
% hold on
% scatter(T2Metric, R2starMedian, [], ADCvec, 'filled');
% %plot(T2Metric, R2starMedian, 'bo')
% ylabel('R2* (s^{-1})')
% xlabel('T2 metric')
% lgd = legend('Tumor', 'Location', 'NorthEast');
% title(lgd, str)
% c = colorbar;
% c.Label.String = 'ADC (mm^2/s)';
% shading interp;
% jet(256);
% set(gca, 'FontSize', 11)
% 
% figure()
% hold on
% scatter(T2Metric, R2starMedian, [], HF_DWIvec, 'filled');
% %plot(T2Metric, R2starMedian, 'bo')
% ylabel('R2* (s^{-1})')
% xlabel('T2 metric')
% lgd = legend('Tumor', 'Location', 'NorthEast');
% title(lgd, str)
% c = colorbar;
% c.Label.String = 'HF_DWI';
% shading interp;
% jet(256);
% set(gca, 'FontSize', 11)
% 
% %%%%%%% 
% [r, pval] = corr(ADCvec, T2Metric, 'type', 'Spearman');
% str = strcat('n=', num2str(length(T2Metric)), ...
%     ', r=', num2str(round(r,3)), ...
%     ', P=', num2str(round(pval,6)));
% figure()
% hold on
% scatter(T2Metric, ADCvec, [], HF_DWIvec, 'filled');
% ylabel('ADC (mm^2/s)')
% xlabel('T2 metric')
% lgd = legend('Tumor', 'Location', 'NorthEast');
% title(lgd, str)
% c = colorbar;
% c.Label.String = 'HF_DWI';
% shading interp;
% jet(256);
% set(gca, 'FontSize', 11)
% 
% [r, pval] = corr(fBVvec, T2Metric, 'type', 'Spearman');
% str = strcat('n=', num2str(length(T2Metric)), ...
%     ', r=', num2str(round(r,3)), ...
%     ', P=', num2str(round(pval,6)));
% figure()
% hold on
% scatter(T2Metric, fBVvec, [], HF_DWIvec, 'filled');
% ylabel('fBV')
% xlabel('T2 metric')
% lgd = legend('Tumor', 'Location', 'NorthEast');
% title(lgd, str)
% c = colorbar;
% c.Label.String = 'HF_DWI';
% shading interp;
% jet(256);
% set(gca, 'FontSize', 11)
% 
% [r, pval] = corr(HF_DWIvec, T2Metric, 'type', 'Spearman');
% str = strcat('n=', num2str(length(T2Metric)), ...
%     ', r=', num2str(round(r,3)), ...
%     ', P=', num2str(round(pval,6)));
% figure()
% hold on
% plot(T2Metric, HF_DWIvec, 'bo')
% ylabel('HF_DWI')
% xlabel('T2 metric')
% lgd = legend('Tumor', 'Location', 'NorthEast');
% title(lgd, str)
% 
% 
% 
% 
% 

clear folderList folderpath foldername i str lgd

end