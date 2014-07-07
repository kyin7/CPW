function plotresult(Psi, Psi_hat, J_Psi, L, w, mu, alpha)
[N p] = size(Psi);
for k=1:p
    psi_hat = Psi_hat(:,k);
    psi = Psi(:,k);
    J = J_Psi(:,k);
    plot_result(psi, psi_hat, J, L, w, mu, alpha, N, k);
end