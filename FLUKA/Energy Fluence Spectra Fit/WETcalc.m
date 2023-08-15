clear all;
close all;
fclose('all');
clc;

%% WET

% 15.22 MeV proton, 103.1 cm distance
% {[Pos 1], [Pos 2], [Pos 3], [Pos 4], [Pos 5]}
% {[IC para], [IC LID], [IC LID], [IC LID], [IC LID]}
E_struct       = {[8.48 8.89],           [1.40 6.81],           [0.97 6.63],           [0.73 6.50],           [0.58 6.41]};
SP_m_struct    = {[5.631E+01 5.468E+01], [2.272E+02 6.124E+01], [2.956E+02 6.256E+01], [3.604E+02 6.356E+01], [4.214E+02 6.426E+01]};
SP_w_struct    = {[5.212E+01 5.018E+01], [2.056E+02 6.205E+01], [2.664E+02 6.338E+01], [3.233E+02 6.438E+01], [3.761E+02 6.510E+01]};
rho_m_struct   = {[9.40E-01 9.30E-01],   [9.40E-01 1.06E+00],   [9.40E-01 1.06E+00],   [9.40E-01 1.06E+00],   [9.40E-01 1.06E+00]}; % in g/cm3
t_m_struct     = {[0.003 0.0130],        [0.003 0.1050],        [0.003 0.1050],        [0.003 0.1050],        [0.003 0.1050]}; % in cm

for i = 1:length(E_struct)
    
    E = E_struct{i}; 
    SP_m = SP_m_struct{i};
    SP_w = SP_w_struct{i};
    rho_m = rho_m_struct{i};
    t_m = t_m_struct{i};
    
    WET = 0;
    
    for j = 1:length(E)
    
        t_w = t_m(j) * rho_m(j) * (SP_m(j)/SP_w(j));
        WET = WET + t_w;
        
    end
    WET
    
end