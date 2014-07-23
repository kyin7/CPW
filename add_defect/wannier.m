function [Psi, Psi_hat, J_Psi, Psi_energy] = wannier(opt, max_level)
N = opt.N;
max_iter = opt.max_iter;
L = opt.L;
Psi = zeros(N,max_level);
Psi_hat = zeros(N,max_level);
J_Psi = zeros(max_iter,max_level);
Psi_energy = zeros(max_level,1);
for level=1:max_level
    tic
%     [psi, psi_hat, J_psi] = wannier_L1(level, opt, Psi_hat);
    [psi, psi_hat, J_psi, J_min] = wannier_L1_randn(level, opt, Psi_hat);
    toc
    Psi(:,level) = psi;
    Psi_hat(:,level) = psi_hat;
    J_Psi(:,level) = J_psi;
    Psi_energy(level,1) = J_min;
end
Psi = sqrt(N/L)*Psi;
Psi_hat = sqrt(N/L)*Psi_hat;

