clear all;
close all;
clc;

baseName = 'DLTS_A_2E10_';
legendName = cell(6,1);
lineSpec = {'b-','r-','g-','k-', 'm-', 'c-'};
counter = 1;
% Plot the DLTS-spectrum for i min into the annealing process
for i = [0 5 10 20 40 80] % annealing time [min]
    fileName = [baseName, num2str(i), 'min.txt'];
    data = load(fileName);
    T = data(:,1);
    C_r = data(:,4);
    w6 = data(:,11);
    hold on
    plot(T, w6./(0.102.*C_r), lineSpec{counter});
    legendName{counter} = [num2str(i), ' min'];    
    counter = counter + 1;
end
xlabel(['Temperature [K]'],'interpreter','latex','FontSize',13)
ylabel(['$DLTS/0.102C_r\,$'],'interpreter','latex','FontSize',13)
legend(legendName)

% Values of (S_max/0.102*C_r) for every defect (peak)
% The vectors are arranged by [0 min; 5 min; ...; 80 min]
V2O_max_1st = [0.01253; 0.01087; 0.009647; 0.007903; 0.005444; 0.002362];
V2O_max_2nd = [0.01321; 0.0124; 0.01089; 0.009691; 0.007505];
new_defect = [0.004161; 0.004296; 0.00463; 0.004688; 0.005372; 0.006817];
t_annealing = [0; 5; 10; 20; 40;  80]; % annealing time [min]
N_d = 3.0*10^13; % doping concentration [cm^-3]

for defect_consentration = [V2O_max_1st V2O_max_2nd new_defect]
    figure()
    N_T = N_d.*defect_consentration;
    plot(t_annealing, N_T, '*')
    p = polyfit(t_annealing, N_T, 1);
    yfit = polyval(p, t_annealing);
    hold on
    plot(t_annealing, yfit, '-ro')
    xlabel(['Time $[min]$'],'interpreter','latex','FontSize',13)
    ylabel(['$N_T\,[cm^{-3}]$'],'interpreter','latex','FontSize',13)
    str = sprintf('Linefit: %0.5g + %0.5g t', p(2), p(1));
    legend('Data', str);
end