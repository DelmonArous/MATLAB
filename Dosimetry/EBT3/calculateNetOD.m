function [netOD] = calculateNetOD(PV_before, PV_after, PV_bckg)

netOD = max(0, log10( (PV_before - PV_bckg) ./ (PV_after - PV_bckg) ));

end