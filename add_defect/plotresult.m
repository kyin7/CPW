function plotresult(Psi, Psi_hat, J_Psi, opt)
[p] = size(Psi,2);
for k=1:p
    psi_hat = Psi_hat(:,k);
    psi = Psi(:,k);
    J = J_Psi(:,k);
    plot_result(psi, psi_hat, J, opt, k);
end