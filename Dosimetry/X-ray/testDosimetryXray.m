clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')


current_mA = [1 1 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 8 8 9 9 10 10 11 11 12 12];
MC_nC = [0.04 0.04 0.10 0.09 0.14 0.14 0.12 0.19 0.19 0.18 0.24 0.24 0.22 ...
    0.29 0.28 0.28 0.33 0.33 0.39 0.39 0.44 0.44 0.50 0.49 0.55 0.55 ...
    0.60 0.60];
TC_nC = [166.4 166.5 460.8 460.7 754.4 755.6 757.0 1051 1051 1051 1338 ...
    1340 1335 1491 1488 1492 1487 1482 1521 1523 1515 1516 1523 1520 ...
    1516 1517 1524 1516];
ratio = MC_nC ./ TC_nC;

figure();
hold on
yyaxis left
plot(current_mA, (MC_nC), 'o', current_mA, (TC_nC), 'x')
ylabel('Chamber reading (nC)')
yyaxis right
plot(current_mA, ratio, 'd')
ylabel('MC/TC')
set(gca, 'FontSize', 14)
hold off
xlabel('Current (mA)')
xlim([min(current_mA)-0.5 max(current_mA)+0.5])
title(legend('MC', 'TC', 'Ratio'), '220 kV, 2 min')

irrad_time_min = [1.0 1.0 1.5 1.5 2.0 2.0 2.5 2.5 3.0 3.0];
MC_nC = [0.23 0.23 0.36 0.36 0.50 0.49 0.61 0.61 0.74 0.75];
TC_nC = [764.1 760.3 1139 1136 1523 1520 1892 1892 2270 2280];
ratio = MC_nC ./ TC_nC;

figure();
hold on
yyaxis left
plot(irrad_time_min, (MC_nC), 'o', irrad_time_min, (TC_nC), 'x')
ylabel('Chamber reading (nC)')
yyaxis right
plot(irrad_time_min, ratio, 'd')
ylabel('MC/TC')
set(gca, 'FontSize', 14)
hold off
xlabel('Irradiation time (min)')
xlim([min(irrad_time_min)-0.5 max(irrad_time_min)+0.5])
title(legend('MC', 'TC', 'Ratio'), '220 kV, 10 mA')

