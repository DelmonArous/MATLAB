function [] = R2starPixelAnalysis(R2starpatients)

r_array_R2star_fBV = zeros(length(R2starpatients), 1);
pval_array_R2star_fBV = zeros(length(R2starpatients), 1);
r_array_R2star_ADC = zeros(length(R2starpatients), 1);
pval_array_R2star_ADC = zeros(length(R2starpatients), 1);

for i = 1:length(R2starpatients)
    
    temp_R2starvaluesvec = R2starpatients{i}.ROIR2starvalues;
    temp_fBVvaluesvec = R2starpatients{i}.ROIfBVvalues;
    ind = find(temp_fBVvaluesvec < 10^(-5)).';
    temp_R2starvaluesvec(ind) = [];
    temp_fBVvaluesvec(ind) = [];
    
%     h = figure();
%     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%     set(h, 'defaultAxesColorOrder', [[0 0 1]; [1 0 0]]);
%     hold on
%     yyaxis left
%     scatter(temp_R2starvaluesvec, temp_fBVvaluesvec, ...
%         15, 'blue', 'filled', 'MarkerEdgeColor', 'black')
%     xlabel('R2* (s^{-1})')
%     ylabel('fBV')
    
    % Correlation
    [r_R2star_fBV, pval_R2star_fBV] = ...
        corr(temp_R2starvaluesvec.', temp_fBVvaluesvec.', 'type', 'Pearson');
    r_array_R2star_fBV(i) = r_R2star_fBV; 
    pval_array_R2star_fBV(i) = pval_R2star_fBV;
    
%     yyaxis right
%     scatter(R2starpatients{i}.ROIR2starvalues, R2starpatients{i}.ROIADCvalues, ...
%         15, 'red', 'filled', 'MarkerEdgeColor', 'black')
%     ylabel('ADC (mm^2/s)')
%     hold off
%     
%     set(gca, 'FontSize', 16)
%     str = {R2starpatients{i}.Patientname, ...
%         strcat('HF_{DWI}=', num2str(round(R2starpatients{i}.HF_DWImedian,3))), ...
%         strcat('HS_{Pimo}=', num2str(round(R2starpatients{i}.HS_Pimomedian,3)))};
%     str{1} = ['\bf ', str{1}, ' \rm'];
%     t = annotation('textbox', 'String', str, 'FitBoxToText', 'on');
%     t.Position = [0.75 0.8 0.1 0.1];
%     t.FontSize = 14;
%     saveas(h, sprintf('pixelplot_R2starVsfBVVsADC_%s.png', ...
%         R2starpatients{i}.Patientname));
    
    % Correlation
    [r_R2star_ADC, pval_R2star_ADC] = ...
        corr(R2starpatients{i}.ROIR2starvalues.', R2starpatients{i}.ROIADCvalues.', 'type', 'Pearson');
    r_array_R2star_ADC(i) = r_R2star_ADC;
    pval_array_R2star_ADC(i) = pval_R2star_ADC;
    
end

%% Histogram plot
figure();
h = histogram(r_array_R2star_fBV);
h.NumBins  = 15;
hold on
h = histogram(r_array_R2star_ADC);
h.NumBins  = 15;
xlim([-1 1])
xlabel('Pearson´s r')
ylabel('No. of patients')
title('Pixel-wise correlation for R_2^*')
legend('R_2^* vs. fBV', 'R_2^* vs. ADC')
set(gca, 'FontSize', 14)
hold off

%% Compare hypoxia status 

%%% Cohort 1
%%%%% R2starpatients{13} = Patient17 (Less hypoxic)
%%%%% R2starpatients{1} = Patient2 (More hypoxic)

%%%%% R2starpatients{36} = Patient42 (Less hypoxic)
%%%%% R2starpatients{43} = Patient49 (More hypoxic)

%%% Cohort 2
%%%%% R2starpatients{70} = Patient84 (Less hypoxic)
%%%%% R2starpatients{91} = Patient111 (More hypoxic)

%%%%% R2starpatients{83} = Patient99 (Less hypoxic)
%%%%% R2starpatients{95} = Patient118 (More hypoxic)

hyp_less = [13 36 70 83];
hyp_more = [1 43 91 95];

% for i = 1:length(hyp_less)
%     
%     h = figure();
%     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%     hold on
%     
%     temp_ADCvaluesvec = R2starpatients{hyp_less(i)}.ROIADCvalues;
%     temp_fBVvaluesvec = R2starpatients{hyp_less(i)}.ROIfBVvalues;
%     ind = find(temp_fBVvaluesvec < 10^(-5)).';
%     temp_ADCvaluesvec(ind) = [];
%     temp_fBVvaluesvec(ind) = [];
%     scatter(temp_ADCvaluesvec, temp_fBVvaluesvec, ...
%         25, 'blue', 'filled', 'MarkerEdgeColor', 'black')
%     
%     temp_ADCvaluesvec = R2starpatients{hyp_more(i)}.ROIADCvalues;
%     temp_fBVvaluesvec = R2starpatients{hyp_more(i)}.ROIfBVvalues;
%     ind = find(temp_fBVvaluesvec < 10^(-5)).';
%     temp_ADCvaluesvec(ind) = [];
%     temp_fBVvaluesvec(ind) = [];
%     scatter(temp_ADCvaluesvec, temp_fBVvaluesvec, ...
%         25, 'red', 'filled', 'MarkerEdgeColor', 'black')
%     
%     xlabel('ADC (mm^2/s)')
%     ylabel('fBV')
%     legend({[R2starpatients{hyp_less(i)}.Patientname ' (Less hypoxic)'], ...
%         [R2starpatients{hyp_more(i)}.Patientname ' (More hypoxic)']});
%     set(gca, 'FontSize', 16)
%     hold off
%     
%     saveas(h, sprintf('ADCfBV_%s_%s.png', ...
%         R2starpatients{hyp_less(i)}.Patientname, R2starpatients{hyp_more(i)}.Patientname));
% 
% end

end