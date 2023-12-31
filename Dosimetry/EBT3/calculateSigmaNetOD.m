function [sigma_netOD] = calculateSigmaNetOD(PV_before, PV_after, ...
    PV_bckg, sigma_PV_before, sigma_PV_after, sigma_bckg)

% sigma_netOD = (1/log(10)) .* sqrt( ...
%     ((sigma_PV_before.^2)./(PV_before - PV_bckg).^2) + ...
%     ((sigma_PV_after.^2)./(PV_after - PV_bckg).^2) + (sigma_bckg.^2) .* ...
%     ((PV_before - PV_after) ./ ((PV_before - PV_bckg) .* ...
%     (PV_after - PV_bckg))).^2);

sigma_netOD = (1/log(10)) .* sqrt( ...
    ((sigma_PV_before.^2)./(PV_before - PV_bckg).^2) + ...
    ((sigma_PV_after.^2)./(PV_after - PV_bckg).^2) + (sigma_bckg.^2) .* ...
    ((PV_before - PV_after) ./ ((PV_before - PV_bckg) .* ...
    (PV_after - PV_bckg))).^2);

end