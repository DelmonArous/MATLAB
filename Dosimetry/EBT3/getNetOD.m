function [netOD, sigma_netOD] = getNetOD(EBT3, calib)

netOD = calculateNetOD(calib.PV_control, calib.PV_bckg, EBT3.PV);

sigma_netOD = calculateSigmaNetOD();


netOD = log10( (calib.PV_control - calib.PV_bckg) / ...
    (EBT3.PV - calib.PV_bckg) );

sigma_netOD = (1/log(10)) * ...
    sqrt( calib.sigma_control^2/(calib.PV_control-calib.PV_bckg)^2 + ...
    EBT3.sigma^2/(EBT3.PV-calib.PV_bckg)^2 + calib.sigma_bckg^2 * ...
    ((calib.PV_control - EBT3.PV) / ...
    ((calib.PV_control-calib.PV_bckg)*(EBT3.PV - calib.PV_bckg)))^2 ); 

end