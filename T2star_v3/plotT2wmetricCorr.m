function [] = plotT2wmetricCorr(R2starpatients)

%% 
parameter = {'R2star', 'ADC', 'fBV', 'HF_DWI', 'HS_Pimo', 'BVD', 'CD'};
labelspec = {'R2* (s^{-1})', 'ADC (mm^2/s)', 'fBV', 'HF_{DWI}', 'HS_{Pimo}', 'BVD', 'CD'};

for i = 1:length(parameter)
    
    valuesvec = []; T2wmetricvec = []; HS_Pimovec = [];
    if strcmp(parameter{i}, 'HS_Pimo')
        for j = 1:length(R2starpatients)
            valuesvec = [valuesvec; ...
                R2starpatients{j}.(sprintf('%smedian', parameter{i}))];
            T2wmetricvec = [T2wmetricvec; R2starpatients{j}.T2wmetric];
            
%             R2starpatients{j}.Patientname
%             R2starpatients{j}.slicelocation
%             R2starpatients{j}.SIfat
%             [R2starpatients{j}.R2starmedian R2starpatients{j}.T2wmetric]
            
        end
%         temp_valuesvec = valuesvec;
%         temp_T2wmetricvec = T2wmetricvec;
%         ind1 = find(isnan(temp_valuesvec)).';
%         temp_valuesvec([ind1]) = [];
%         temp_T2wmetricvec([ind1]) = [];
%         [r, pval] = corr(temp_T2wmetricvec, temp_valuesvec, 'type', 'Spearman'); % HF_DWIvec
%         str = strcat('n=', num2str(length(temp_valuesvec)), ...
%             ', R^2=', num2str(round(r^2,3)), ...
%             ', P=', num2str(round(pval,3)));
        
        ind1 = find(isnan(valuesvec)).';
        ind2 = find(isnan(T2wmetricvec)).';
        valuesvec([ind1, ind2]) = [];
        T2wmetricvec([ind1, ind2]) = [];
    else
        for j = 1:length(R2starpatients)
            valuesvec = [valuesvec; ...
                R2starpatients{j}.(sprintf('%smedian', parameter{i}))];
            T2wmetricvec = [T2wmetricvec; R2starpatients{j}.T2wmetric];
            HS_Pimovec = [HS_Pimovec; R2starpatients{j}.HS_Pimomedian];
        end
%         temp_valuesvec = valuesvec;
%         temp_T2wmetricvec = T2wmetricvec;
%         ind1 = find(isnan(temp_valuesvec)).';
%         temp_valuesvec([ind1]) = [];
%         temp_T2wmetricvec([ind1]) = [];
        
        ind1 = find(isnan(valuesvec)).';
        ind2 = find(isnan(T2wmetricvec)).';
        ind3 = find(isnan(HS_Pimovec)).';
        valuesvec([ind1, ind2, ind3]) = [];
        T2wmetricvec([ind1, ind2, ind3]) = [];
        HS_Pimovec([ind1, ind2, ind3]) = [];
    end
        
    [r, pval] = corr(T2wmetricvec, valuesvec, 'type', 'Spearman'); % HF_DWIvec
    X = [ones(length(T2wmetricvec), 1) T2wmetricvec];
    beta = X\valuesvec;
    str = strcat('n=', num2str(length(valuesvec)), ...
        ', r=', num2str(round(r,3)), ...
        ', P=', num2str(round(pval,9)));
    
    %% Plot
    h = figure();
    hold on
    
    if strcmp(parameter{i}, 'HS_Pimo')
        %%% T2metric against HS_Pimo without colorbar
        plot(T2wmetricvec, valuesvec, 'ko') % , T2wmetricvec, X*beta, 'r-');
    else
        %%% R2* against ADC, fBV and HS_Pimo with colorbar
        scatter(T2wmetricvec, valuesvec, [], HS_Pimovec, 'filled', ...
            'MarkerEdgeColor', 'black');
        %plot(T2wmetricvec, X*beta, 'r-');
        c = colorbar;
        c.Label.String = 'HS_{Pimo}';
        shading interp;
        jet(256);
    end
    
    set(gca, 'FontSize', 14)
    xlabel('T2 metric')
    ylabel(labelspec{i})
    lgd = legend('Median tumor value', 'Location', 'best'); 
    %%% num2str(10^(floor(log(abs(pval_Q2))./log(10)) + 1)));
    title(lgd, str)
    %saveas(h, strcat('T2metricCorr_', parameter{i}, '.png'))
    
end

end