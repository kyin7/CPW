function [Psi, Psi_hat, J_Psi] = cpws(L, w, mu, alpha, lambda, r,  N, max_level, max_iter)
Psi = zeros(N,max_level);
Psi_hat = zeros(N,max_level);
J_Psi = zeros(max_iter,max_level);
h = floor(w/L*N);
mu_s = mu * (2*pi/L)^2*sqrt(N/L);
lambda_s = lambda/(2*pi/L)^2;
r_s = r/(2*pi/L)^2;
for level=1:max_level
    tic
    [psi, psi_hat, J_psi] = cpw(lambda_s*level^alpha, r_s*level^alpha, mu_s/level^alpha, N, h, level, Psi_hat, max_iter);
    toc
    Psi(:,level) = psi;
    Psi_hat(:,level) = psi_hat;
    J_Psi(:,level) = J_psi;
end
Psi = sqrt(N/L)*Psi;
Psi_hat = sqrt(N/L)*Psi_hat;
