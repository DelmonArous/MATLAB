function [R2starpatients] = R2star_fBV_hypoxiaRegistration(R2starpatients)

%% Cohort 1 + Cohort 2
xvec = [1 1 24 0 20 5 19 -8 -1 4 -8 17 -3 23 -3 ...
    -5 10 9 0 4 4 22 -6 -19 ...
    23 -1 10 13 24 8 14 15 11 ...
    26 13 -5 9 6 10 10 -6 13 ...
    4 2 25 11 ...
    1 14 19 2 12 8 4 8 ... % Cohort 2
    -9 9 12 29 10 16 17 3 27 ...
    1 2 11 25 1 2 5 2 4 ...
    2 24 3 4 6 12 -6 28 12 ...
    9 10 10 6 7 6 2 4 33 4 ...
    13 -11 -7 -11 11 20];
yvec = [1 1 -23 1 -18 -4 -18 9 2 -3 10 -17 5 -21 5 ...
    6 -10 -7 2 -3 -2 -20 7 21 ...
    -22 2 -8 -10 -22 -6 -14 -13 -10 ...
    -24 -12 6 -7 -6 -8 -9 7 -12 ...
    -3 -1 -25 -10 ...
    1 -13 -18 -2 -10 -7 -2 -6 ... % Cohort 2
    9 -9 -10 -28 -8 -16 -16 -2 -27 ...
    0 -1 -10 -24 -1 -1 -4 -2 -2 ...
    0 -23 -2 -4 -4 -11 7 -26 -11 ...
    -8 -10 -9 -4 -6 -4 -2 -2 -33 -2 ...
    -11 12 8 12 -11 -18];

% xvec = [1]; yvec = [1];

%%
R2star_thres = 5:0.5:34;
fBV_thres = 0.01:0.01:0.2;
cc = jet(length(fBV_thres));
ColorSet = varycolor(length(fBV_thres));

for i = 1:length(R2starpatients)
    
    %% Resize fBV image to equal R2* image size
    tempmask = zeros(size(R2starpatients{i}.R2starMap));
    s = size(R2starpatients{i}.fBVmap);
    
    x = abs(R2starpatients{i}.ImagePositionPatient(1) - ...
        R2starpatients{i}.R2starImagePositionPatient(1)) / ...
        R2starpatients{i}.PixelSpacing(1);
    y = abs(R2starpatients{i}.ImagePositionPatient(2) - ...
        R2starpatients{i}.R2starImagePositionPatient(2)) / ...
        R2starpatients{i}.PixelSpacing(2);
    
    x = x + xvec(i); y = y + yvec(i);
    tempmask(x:x+s(1)-1, y:y+s(2)-1) = R2starpatients{i}.fBVmap;
    ROIfBVmap = (tempmask + 1e-6) .* R2starpatients{i}.R2starROImask;
    
    %% Define hypoxic region through fBV threshold at 0.03 (Cohort 1) or 0.17 (Cohort 2), or 0.02 (Cohort 1+2)
    for k = 1:length(fBV_thres)
    
        ROIfBVmask_hypoxic = ROIfBVmap > 0 & ROIfBVmap <= fBV_thres(k);
        ROIfVBmap_hypoxic = ROIfBVmap .* ROIfBVmask_hypoxic;
        R2starpatients{i}.HF_fBV(k) = nnz(ROIfVBmap_hypoxic) / nnz(ROIfBVmap);
        
        %% Define residual hypoxic region through R2* threshold
        resid_mask = (ROIfBVmask_hypoxic == 0);
        resid_R2starMap = R2starpatients{i}.R2starMap .* resid_mask;
        resid_ROIR2starMap = resid_R2starMap .* R2starpatients{i}.R2starROImask;
        
        %     h = figure();
        %     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %     imshow(ROIfBVmap, [min(ROIfBVmap(:)) max(ROIfBVmap(:))], 'colormap', jet(256))
        %     title([R2starpatients{i}.Patientname ' fBV ROI'])
        %
        %     h = figure();
        %     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %     imshow(ROIfVBmap_hypoxic, [min(ROIfBVmap(:)) max(ROIfBVmap(:))], 'colormap', jet(256))
        %     title([R2starpatients{i}.Patientname ' fBV ROI hypoxia'])
        %
        %     h = figure();
        %     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %     imshow(resid_ROIR2starMap, [0 50], 'colormap', jet(256))
        %     title([R2starpatients{i}.Patientname ' residual R2* ROI hypoxia'])
        
        %% Loop through every R2* threshold defining hypoxia
        resid_ROIR2starvalues = reshape(resid_ROIR2starMap, 1, []); % resid_ROIR2starMap(:);
        resid_ROIR2starvalues(resid_ROIR2starvalues == 0) = [];
        resid_ROIR2starvalues(isnan(resid_ROIR2starvalues)) = [];
        resid_ROIR2starvalues(isinf(resid_ROIR2starvalues)) = [];
        
        for j = 1:length(R2star_thres)
            
            resid_ROIR2starvalues_hypoxic = resid_ROIR2starvalues(resid_ROIR2starvalues >= R2star_thres(j));
            R2starpatients{i}.HF_R2star(j) = length(resid_ROIR2starvalues_hypoxic) / ...
                nnz(R2starpatients{i}.R2starROImask);
            R2starpatients{i}.HF_tot(j,k) = R2starpatients{i}.HF_R2star(j) + ...
                R2starpatients{i}.HF_fBV(k);
            
        end
    
    end
     
%     R2starpatients{i}.Patientname
%     length(R2starpatients{i}.ROIR2starvalues)
%     nnz(ROIfBVmap)
%     nnz(R2starpatients{i}.R2starROImask)
  
end
 
% %%
r_array = zeros(length(R2star_thres), length(fBV_thres));
pval_array = zeros(length(R2star_thres), length(fBV_thres));

for k = 1:length(fBV_thres)
    
    r_array = [];
    pval_array = [];
    
    for j = 1:length(R2star_thres)
        
        HF_totvec = []; HF_DWIvec = []; HS_Pimovec = [];
        
        for i = 1:length(R2starpatients)
            
            HF_totvec = [HF_totvec; R2starpatients{i}.HF_tot(j,k)];
            HS_Pimovec = [HS_Pimovec; R2starpatients{i}.HS_Pimomedian];
            
        end
        
        ind1 = find(isnan(HF_totvec)).';
        ind2 = find(isnan(HS_Pimovec)).';
        HF_totvec([ind1, ind2]) = [];
        HS_Pimovec([ind1, ind2]) = [];
        
        [r, pval] = corr(HF_totvec, HS_Pimovec, 'type', 'Pearson');
        % X = [ones(length(HS_Pimovec), 1) HS_Pimovec];
        % beta = X\HF_totvec; % HF_DWIvec
        r_array(j,k) = r; pval_array(j,k) = pval;
        
    end
    
    h = figure(1);
    hold on
    plot(R2star_thres, r_array(:,k), 'Color', cc(k,:))
    xlabel('R2*_{thres} (s^{-1})')
    ylabel('Pearson´s r')
    xlim([R2star_thres(1) R2star_thres(end)])
    xticks(R2star_thres(1):5:R2star_thres(end))
    c = colorbar;
    c.Label.String = 'fVB_{thres} (a.u)';
    %str = strcat('fBV_{thres}=', num2str(fBV_thres(k)));
    
%     legend(str, 'Location', 'NorthEastOutside')
%     set(h, 'Colormap', ColorSet);
%     c = colorbar;
%     c.Label.String = 'fVB_{thres} (a.u)';
%     set(c, 'ylim', [fBV_thres(1) fBV_thres(end)])
    %set(c, 'YDir', 'reverse');
    set(gca, 'FontSize', 14)     
    hold off
    
    h = figure(2);
    hold on
    plot(R2star_thres, pval_array(:,k), 'Color', cc(k,:))
    plot(R2star_thres, repmat(0.05, 1, length(R2star_thres)), 'k--', 'LineWidth', 1)
    xlabel('R2*_{thres} (s^{-1})')
    ylabel('P-value (Pearson)')
    xlim([R2star_thres(1) R2star_thres(end)])
    xticks(R2star_thres(1):5:R2star_thres(end))
    %legend(str, 'Location', 'NorthEastOutside')
    c = colorbar;
    c.Label.String = 'fVB_{thres} (a.u)';
    
%     set(h, 'Colormap', ColorSet);
%     c = colorbar;
%     c.Label.String = 'fVB_{thres} (a.u)';
%     set(c, 'ylim', [fBV_thres(1) fBV_thres(end)])
%     c = colorbar;
%     c.Label.String = 'fVB_{thres} (a.u)';
%     set(c, 'ylim', [fBV_thres(1) fBV_thres(end)])
%     set(c, 'YDir', 'reverse');
    set(gca, 'FontSize', 14)
    hold off
    
    %% Plot
%     h = figure();
%     hold on
%     
%     %%%
%     %plot(HS_Pimovec, HF_totvec, 'ko', HS_Pimovec, X*beta, 'r-')
%     scatter(HS_Pimovec, HF_totvec, [], HF_DWIvec, 'filled', ...
%         'MarkerEdgeColor', 'black');
%     %plot(HS_Pimovec, X*beta, 'r-');
%     c = colorbar;
%     c.Label.String = 'HF_{DWI}';
%     shading interp;
%     jet(256);
%     
%     xlabel('HS_{Pimo}')
%     ylabel('HF_{R2*}+HF_{fBV}')
%     str = ['fBV_{thres} = 0.03, R2*_{thres} = ', num2str(thres(j)), ' s^{-1}'];
%     lgd = legend(str, 'Location', 'SouthEast');
%     str = strcat('n=', num2str(length(HF_totvec)), ...
%         ', r=', num2str(round(r,3)), ...
%         ', P=', num2str(round(pval,3)));
%     title(lgd, str)
%     set(gca, 'FontSize', 14)
%     hold off
    
end

% lgd = cell(length(fBV_thres),1);
% for i = 1:length(fBV_thres)
%     lgd{i} = strcat('fBV_{thres}=', num2str(fBV_thres(i)));
% end
% leg = legend(lgd, 'Location', 'NorthEastOutside');
% set(leg, 'FontSize', 11)

%%
% h = figure();
% 
% yyaxis left
% plot(R2star_thres, r_array(:,1), '-o')
% xlabel('R2*_{thres} (s^{-1})')
% ylabel('Pearson´s r')
% xlim([R2star_thres(1) R2star_thres(end)])
% xticks(R2star_thres(1):5:R2star_thres(end))
% 
% yyaxis right
% plot(R2star_thres, pval_array(:,1), '-d', ...
%     R2star_thres, repmat(0.05, 1, length(R2star_thres)), '--')
% ylabel('P-value')
% 
% title('HF_{R2*}+HF_{fBV}-HS_{Pimo} correlation')
% set(gca, 'FontSize', 14)
% saveas(h, 'HF_R2starPlusHF_fBV_HS_Pimo_PearsonCorr.png')

end