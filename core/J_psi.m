function J = J_psi(mu, N, psi, psi_hat)
% J(psi) = <psi, H psi> + 1/mu*|psi|_1
J = sum(0.5/N*((-N/2:N/2-1)'.^2).*(abs(psi_hat).^2) + abs(psi)/mu);