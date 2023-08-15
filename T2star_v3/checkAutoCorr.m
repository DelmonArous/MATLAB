function [] = checkAutoCorr(R2starpatients)

% R2starvec = [];
% prctilevec = [];
% for j = 1:99
%     R2starvec = [R2starvec; R2starpatients{1}.R2starpercentile(j)];
%     prctilevec = [prctilevec; j];
% end
% [r, lags] = xcov(R2starvec, 'coeff');
% mar = max(r(lags > 0));
% figure()
% hold on
% plot(lags, r, 'b.-')
% xlim([0 100])
% xlabel('Lags')
% ylabel('Correlation')

for i = 2:length(R2starpatients)
    
    R2starvec = [];
    prctilevec = [];
    
    for j = 1:99
        R2starvec = [R2starvec; R2starpatients{i}.R2starpercentile(j)];
        prctilevec = [prctilevec; j];
    end
    
    ind1 = find(isnan(R2starvec)).';
    ind2 = find(isinf(R2starvec)).';
    R2starvec([ind1, ind2]) = [];
    prctilevec([ind1, ind2]) = [];
    
%     [rb, lags] = xcov(R2starvec, 'coeff');
%     marb(i-1) = max(rb(lags > 0));
    
    
    figure();
    autocorr(R2starvec, length(R2starvec)-1)
    xlabel('Lags')
    ylabel('Correlation')
    title('Autocorrelation Function with 95% Confidence Interval')
    
%     [h, p, Qstat, crit] = lbqtest(R2starvec, 'Lags', 0);
%     format long
    % prctilevec(1):10:prctilevec(end)
    %[h; p]
     
end

% pval = mean(marb >= mar)
% rc = quantile(marb, 0.95)
% 
% plot(lags, repmat(rc, 1, length(lags)), 'b--')

end