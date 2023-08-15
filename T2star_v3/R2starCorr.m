function [R2starpatients] = R2starCorr(xlsxpath, R2starpatients)

%% Read xls file
data = readXLSXdocument(xlsxpath);
data_patientnames = data(2:end, 1);
HF_DWI = cell2mat(data(2:end, 4));
HS_Pimo = cell2mat(data(2:end, 5));
BVD = cell2mat(data(2:end, 6));
CD = cell2mat(data(2:end, 7));
N_lymph_status = cell2mat(data(2:end, 8));

%% String variabels
parameter = {'ADC', 'fBV', 'HF_DWI', 'HS_Pimo', 'BVD', 'CD'};
xlabelspec = {'ADC (mm^2/s)', 'fBV', 'HF_{DWI}', 'HS_{Pimo}', 'BVD', 'CD'};
ylabelspec = {'ADC', 'fBV', 'HF_{DWI}', 'HS_{Pimo}', 'BVD', 'CD'};

%%% Cohort 2
% parameter = {'ADC', 'fBV', 'HF_DWI', 'HS_Pimo'};
% xlabelspec = {'ADC (mm^2/s)', 'fBV', 'HF_{DWI}', 'HS_{Pimo}'};
% ylabelspec = {'ADC', 'fBV', 'HF_{DWI}', 'HS_{Pimo}'};

%% Add HF_DWI and HS_pimo score
for i = 1:length(R2starpatients)
    for j = 1:length(data_patientnames)
        
        if strcmp(R2starpatients{i}.Patientname, data_patientnames(j))
            R2starpatients{i}.HF_DWImedian = HF_DWI(j);
            R2starpatients{i}.HS_Pimomedian = HS_Pimo(j);
            R2starpatients{i}.BVDmedian = BVD(j);
            R2starpatients{i}.CDmedian = CD(j);
            R2starpatients{i}.N_lymph_status_PSAmedian = N_lymph_status(j);
        end
        
    end
end

r_array = zeros(99, length(parameter));
pval_array = zeros(99, length(parameter));

for k = 1:99
    
    for i = 1:length(parameter)
        
        valuesvec = []; R2starvec = []; HS_Pimovec = [];
        
        if strcmp(parameter{i}, 'HS_Pimo')
            for j = 1:length(R2starpatients)
                valuesvec = [valuesvec; ...
                    R2starpatients{j}.(sprintf('%smedian', parameter{i}))];
                R2starvec = [R2starvec; ...
                    R2starpatients{j}.R2starpercentile(k)];
            end
            
            ind1 = find(isnan(valuesvec)).';
            ind2 = find(isnan(R2starvec)).';
            ind3 = find(isinf(R2starvec)).';
            valuesvec([ind1, ind2, ind3]) = [];
            R2starvec([ind1, ind2, ind3]) = [];
            
        else
            for j = 1:length(R2starpatients)
                valuesvec = [valuesvec; ...
                    R2starpatients{j}.(sprintf('%smedian', parameter{i}))];
                R2starvec = [R2starvec; ...
                    R2starpatients{j}.R2starpercentile(k)];
                HS_Pimovec = [HS_Pimovec; ...
                    R2starpatients{j}.HS_Pimomedian]; % fBVmedian
            end
            
            ind1 = find(isnan(valuesvec)).';
            ind2 = find(isnan(R2starvec)).';
            ind3 = find(isinf(R2starvec)).';
            %ind4 = find(isnan(HS_Pimovec)).';
            valuesvec([ind1, ind2, ind3]) = [];
            R2starvec([ind1, ind2, ind3]) = [];
            HS_Pimovec([ind1, ind2, ind3]) = [];
        end
        
        %% Linear regression
        [r, pval] = corr(R2starvec, valuesvec, 'type', 'Spearman'); 
        str = strcat('n=', num2str(length(valuesvec)), ...
            ', R^2=', num2str(round(r^2,3)), ...
            ', P=', num2str(round(pval,6)));
        
        X = [ones(length(valuesvec), 1) valuesvec];
        beta = X\R2starvec;
        r_array(k,i) = r; pval_array(k,i) = pval;
        
        %% Plot
%         h = figure();
%         hold on
%         
%         if strcmp(parameter{i}, 'HS_Pimo')
%             % R2* against HS_Pimo without colorbar
%             plot(valuesvec, R2starvec, 'ko', valuesvec, X*beta, 'r-');
%         else
%             % R2* against ADC, fBV and HF_DWI with HS_Pimo as colorbar
%             scatter(valuesvec, R2starvec, [], HS_Pimovec, 'filled', ...
%                 'MarkerEdgeColor', 'black');
%             plot(valuesvec, X*beta, 'r-');
%             c = colorbar;
%             c.Label.String = 'HS_{Pimo}';
%             shading interp;
%             jet(256);
%         end
%         
%         set(gca, 'FontSize', 14)
%         xlabel(xlabelspec{i})
%         ylabel('R2* (s^{-1})')
%         lgd = legend('Median tumor value', 'Location', 'best');
%         %%%%% num2str(10^(floor(log(abs(pval_Q2))./log(10)) + 1)));
%         title(lgd, str)
        %saveas(h, strcat('R2starCorr_', parameter{i}, '.png'))
        
    end
    
end

%%
% pval_array = [];
% for k = 1:99
%     
%     LNneg = []; LNpos = []; 
%     
%     for i = 1:length(R2starpatients)
%         
%         if R2starpatients{i}.N_lymph_status_PSAmedian == 0
%             LNneg = [LNneg; R2starpatients{i}.R2starpercentile(k)];
%         else
%             LNpos = [LNpos; R2starpatients{i}.R2starpercentile(k)];
%         end
%         
%     end
%     
%     pval = ranksum(LNneg, LNpos); 
%     pval_array = [pval_array; pval];
%     
% %     h = figure();
% %     group = [1 * ones(size(LNneg)); 2 * ones(size(LNpos))];
% %     boxplot([LNneg; LNpos], group, 'Colors', 'k', 'Labels', {'LN-', 'LN+'}, ...
% %         'Symbol', 'ko')
% %     set(gca, 'FontSize', 14)
% %     str = strcat('P=', num2str(round(pval, 2)));
% %     t = annotation('textbox', 'String', str, 'FitBoxToText', 'on');
% %     t.Position = [0.70 0.8 0.1 0.1];
% %     t.FontSize = 14;
% %     xlabel('Lymph node status')
% %     ylabel('R2* (s^{-1})')
%      
% end
% 
% h = figure();
% 
% plot(1:99, pval_array, '-rd')
% xlabel('Percentile')
% ylabel('P-value')
% xlim([1 99])
% xticks([1,10,20,30,40,50,60,70,80,90,99])
% title('R2*-Lymph node status')
% set(gca, 'FontSize', 14)


%%
 
for i = 1:size(r_array,2)
    
    h = figure();
    
    hold on
    yyaxis left
    plot(1:99, r_array(:,i), '-o')
    xlabel('Percentile')
    ylabel('Spearman´s r')
    xlim([1 99])
    xticks([1,10,20,30,40,50,60,70,80,90,99])
    %xlim([20 80])
    %xticks([20,30,40,50,60,70,80])
    
    yyaxis right
    plot(1:99, pval_array(:,i), '-d', 1:99, repmat(0.05, 1, 99), '--')
    ylabel('P-value')
    title({'Percentile tumor value', ['R2* vs. ' ylabelspec{i} ' correlation']})
    set(gca, 'FontSize', 14)     
    hold off
    
    saveas(h, strcat('R2starCorrPrctile_Spearman_', parameter{i}, '.png'))
    
end

end