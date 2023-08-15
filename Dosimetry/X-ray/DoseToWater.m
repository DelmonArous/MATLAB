function [D_w, SD_D_w] = DoseToWater(M_u, N_k, K_u, mu_en, P_u, k_TP, SD_N_k) 

D_w     = mean(M_u) * N_k * K_u * mu_en * P_u * k_TP;
SD_D_w  = D_w * sqrt((SD_N_k/N_k)^2 + (std(M_u)/mean(M_u))^2);

end