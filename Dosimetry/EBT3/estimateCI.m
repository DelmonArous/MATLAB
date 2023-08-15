function varargout = estimateCI(x, significance_level)

alpha = [significance_level/2 1-(significance_level/2)];

if size(x,3) == 1
    avg     = mean(x);
    SEM     = std(x) ./ sqrt(length(x));        % Standard Error
    ts      = tinv(alpha, length(x)-1);         % T-Score
    CI      = avg + ts .* SEM;                  % Confidence Intervals
    if nargout >= 1; varargout{1} = avg; end
    if nargout >= 2; varargout{2} = CI; end
elseif size(x,3) > 1
    avg     = mean(x, 3);
    SEM     = std(x, 0, 3) ./ sqrt(size(x,3));  % Standard Error
    ts      = tinv(alpha, size(x,3)-1);         % T-Score
    CI_min  = avg + ts(1) .* SEM;               % Confidence Intervals
    CI_max  = avg + ts(2) .* SEM;
    if nargout >= 1; varargout{1} = avg; end
    if nargout >= 2; varargout{2} = CI_min; end
    if nargout >= 3; varargout{3} = CI_max; end
end

end