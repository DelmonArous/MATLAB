function [b] = LowDoseXrayDosimetry(M_u, N_k, K_u, mu_en, P_u, k_TP, ...
    SD_N_k, t_rad)

% Initialize
D_w_vec = [];

% Estimate dose to water
for i = 1:length(M_u)
    [D_w, SD_D_w] = ...
        DoseToWater(M_u{i}, N_k, K_u, mu_en, P_u, k_TP, SD_N_k);
    D_w_vec = [D_w_vec; D_w];
end

% Linear regression
X = [ones(length(t_rad),1) t_rad];
b = X\D_w_vec;
D_fit = X*b;

% Plot
figure();
plot(t_rad, D_w_vec, 'bo', t_rad, D_fit, '-')
xlabel('Irradiation time (s)')
ylabel('Dose to water D_w (Gy)')
xlim([t_rad(1)-1 t_rad(end)+1])

end