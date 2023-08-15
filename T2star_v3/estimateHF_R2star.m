function [] = estimateHF_R2star(R2starpatients)

R2star_thres = 5:0.5:34;

for i = 1:length(R2starpatients)
    
    ROIR2starvalues = R2starpatients{i}.ROIR2starvalues;
    R2starpatients{i}.Patientname
    
    for j = 1:length(R2star_thres)
        
        temp_ROIR2starvalues = ROIR2starvalues(ROIR2starvalues >= R2star_thres(j));
        temp_HF_R2star = sum(~isnan(temp_ROIR2starvalues) & ~isinf(temp_ROIR2starvalues)) / ...
            sum(~isnan(ROIR2starvalues) & ~isinf(ROIR2starvalues));
        R2starpatients{i}.HF_R2star(j) = temp_HF_R2star;
        
    end
    
end

parameters = {'HF_DWI', 'HS_Pimo'};
labels = {'HF_{DWI}', 'HS_{Pimo}'};

for k = 1:length(parameters)
    
    r_array = []; 
    pval_array = []; 
    DWpval_array = [];
    
    for j = 1:length(R2star_thres)
        
        HF_R2starvec = []; valuesvec = [];
        
        for i = 1:length(R2starpatients)
            
            HF_R2starvec = [HF_R2starvec; R2starpatients{i}.HF_R2star(j)];
            valuesvec = [valuesvec; R2starpatients{i}.(sprintf('%smedian', parameters{k}))];
            
        end
        
        ind1 = find(isnan(HF_R2starvec)).';
        ind2 = find(isnan(valuesvec)).';
        HF_R2starvec([ind1, ind2]) = [];
        valuesvec([ind1, ind2]) = [];
        
        [r, pval] = corr(HF_R2starvec, valuesvec, 'type', 'Pearson');
        X = [ones(length(valuesvec), 1) valuesvec];
        beta = X\HF_R2starvec; % HF_DWIvec
        r_array = [r_array; r]; pval_array = [pval_array; pval];
        
        %% Test for Autocorrelation Among Residuals using a Durbin-Watson test
        [~, ~, resid] = regress(HF_R2starvec, valuesvec);
        [DWpval, ~] = dwtest(resid, valuesvec, 'Exact', 'both');
        DWpval_array = [DWpval_array; DWpval];
        
        %% Plot
%         h = figure();
%         hold on
%         
%         %     %%%
%         %     scatter3(HF_R2starvec, HF_DWIvec, thresh(j), ...
%         %         'filled', 'MarkerEdgeColor', 'black')
%         
%         %%%
%         plot(valuesvec, HF_R2starvec, 'ko') % valuesvec, X*beta, 'r-')
%         
%         xlabel(labels{k})
%         ylabel('HF_{R2*}')
%         str = ['R2*_{thres}=', num2str(R2star_thres(j)), ' s^{-1}'];
%         lgd = legend(str, 'Location', 'SouthEast');
%         str = strcat('n=', num2str(length(HF_R2starvec)), ...
%             ', r=', num2str(round(r,3)), ...
%             ', P=', num2str(round(pval,6)));
%         title(lgd, str)
%         set(gca, 'FontSize', 14)
        
    end
    
    h = figure();
    
    yyaxis left
    plot(R2star_thres, r_array(:,1), '-o')
    xlabel('R2*_{thres} (s^{-1})')
    ylabel('Pearson´s r')
    xlim([R2star_thres(1) R2star_thres(end)])
    xticks(R2star_thres(1):5:R2star_thres(end))
    
    yyaxis right
    plot(R2star_thres, pval_array(:,1), '-d', ...
        R2star_thres, repmat(0.05, 1, length(R2star_thres)), '--')
    ylabel('P-value')
    
    title(['HF_{R2*}-', labels{k}, ' correlation'])
    set(gca, 'FontSize', 14)
    %saveas(h, strcat('HF_R2star_', parameters{k}, '_PearsonCorr.png'))
    
    h = figure();
    plot(R2star_thres, DWpval_array(:,1), '-rd', ...
        R2star_thres, repmat(0.05, 1, length(R2star_thres)), 'r--')
    xlabel('R2*_{thres} (s^{-1})')
    ylabel('P-value')
    xlim([R2star_thres(1) R2star_thres(end)])
    xticks(R2star_thres(1):5:R2star_thres(end))
    
    title({['Test for HF_{R2*}-', labels{k}, ' autocorrelation']; ...
        ' by two-sided Durbin-Watson test'})
    set(gca, 'FontSize', 14)
    
end

%%

fBV_thres = 0.01:0.01:0.7;

for i = 1:length(R2starpatients)
    
    ROIfBVvalues = R2starpatients{i}.ROIfBVvalues;
    
    for j = 1:length(fBV_thres)
        
        temp_ROIfBVvalues = ROIfBVvalues(ROIfBVvalues <= fBV_thres(j));
        temp_HF_fBV = sum(~isnan(temp_ROIfBVvalues) & ~isinf(temp_ROIfBVvalues)) / ...
            sum(~isnan(ROIfBVvalues) & ~isinf(ROIfBVvalues));
        R2starpatients{i}.HF_fBV(j) = temp_HF_fBV;
        
    end
    
end

r_array = [];
pval_array = [];

for j = 1:length(fBV_thres)
    
    HF_fBVvec = []; valuesvec = [];
    
    for i = 1:length(R2starpatients)
        
        HF_fBVvec = [HF_fBVvec; R2starpatients{i}.HF_fBV(j)];
        valuesvec = [valuesvec; R2starpatients{i}.HS_Pimomedian];
        
    end
    
    ind1 = find(isnan(HF_fBVvec)).';
    ind2 = find(isnan(valuesvec)).';
    HF_fBVvec([ind1, ind2]) = [];
    valuesvec([ind1, ind2]) = [];
    
    [r, pval] = corr(HF_fBVvec, valuesvec, 'type', 'Pearson');
    X = [ones(length(valuesvec), 1) valuesvec];
    beta = X\HF_fBVvec; % HF_DWIvec
    r_array = [r_array; r]; pval_array = [pval_array; pval];
    
    %% Plot
%     h = figure();
%     hold on
%     
% % % % % %     %%%
% % % % % %     scatter3(HF_fBVvec, HF_DWIvec, fBV_thres(j), ...
% % % % % %         'filled', 'MarkerEdgeColor', 'black')
%     
%     %%
%     plot(valuesvec, HF_fBVvec, 'ko') %  valuesvec, X*beta, 'r-')
%     
%     xlabel('HS_{Pimo}')
%     ylabel('HF_{fBV}')
%     str = ['fBV_{thres}=', num2str(fBV_thres(j))];
%     lgd = legend(str, 'Location', 'SouthEast');
%     str = strcat('n=', num2str(length(HF_fBVvec)), ...
%         ', r=', num2str(round(r,3)), ...
%         ', P=', num2str(round(pval,6)));
%     title(lgd, str)
%     set(gca, 'FontSize', 14)
    
end

h = figure();

yyaxis left
plot(fBV_thres, r_array(:,1), '-o')
xlabel('fBV_{thres}')
ylabel('Pearson´s r')
xlim([fBV_thres(1) fBV_thres(end)])
xticks([fBV_thres(1) 0.1 0.2 0.3 0.4 0.5 0.6 fBV_thres(end)])

yyaxis right
plot(fBV_thres, pval_array(:,1), '-d', ...
    fBV_thres, repmat(0.05, 1, length(fBV_thres)), '--')
ylabel('P-value')

title('HF_{fBV}-HS_{Pimo} correlation')
set(gca, 'FontSize', 14)
%saveas(h, 'HF_fBV_HS_Pimo_PearsonCorr.png')

end