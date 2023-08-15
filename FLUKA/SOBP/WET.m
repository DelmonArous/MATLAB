clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%%

rho_polyethylene    = 9.40*10^(-1); % in g/cm3
rho_paraffin        = 9.30*10^(-1);
rho_polystyrene     = 1.06*10^(0);
rho_water           = 1.00*10^(0);
rho_air             = 1.20479E-03;

t_IC        = 0.003; % in cm
t_parafilm  = 0.0130;
t_lid       = 0.1050;

%% Pos 1

S_IC_m          = 5.631*10^1; % in MeV cm2/g
S_IC_w          = 5.212*10^1; 

t_IC_w = t_IC * (rho_polyethylene/rho_water) * (S_IC_m / S_IC_w);

S_parafilm_m    = 5.468*10^1;
S_parafilm_w    = 5.018*10^1;

t_parafilm_w = t_parafilm * (rho_paraffin/rho_water) * (S_parafilm_m / S_parafilm_w);

t_IC_w + t_parafilm_w 

%% Pos 2

S_IC_m          = 2.272*10^2;
S_IC_w          = 2.056*10^2;

t_IC_w = t_IC * (rho_polyethylene/rho_water) * (S_IC_m / S_IC_w);

S_lid_m         = 6.124*10^1;
S_lid_w         = 6.205*10^1;

t_lid_w = t_lid * (rho_polystyrene/rho_water) * (S_lid_m / S_lid_w);

t_IC_w + t_lid_w - (7.9 * (rho_air/rho_water) * ((3.712*10^1)/(4.230*10^1)) )


%% Pos 3

S_IC_m          = 2.956*10^2;
S_IC_w          = 2.664*10^2;

t_IC_w = t_IC * (rho_polyethylene/rho_water) * (S_IC_m / S_IC_w);

S_lid_m         = 6.256*10^1;
S_lid_w         = 6.338*10^1;

t_lid_w = t_lid * (rho_polystyrene/rho_water) * (S_lid_m / S_lid_w);

t_IC_w + t_lid_w - (4.6 * (rho_air/rho_water) * ((3.712*10^1)/(4.230*10^1)) )
