function J = J_psi(L, mu, N, psi, psi_hat)
% J(psi) = <psi, H psi> + 1/mu*|psi|_1
% J = sum(0.5/N*((-N/2:N/2-1)'.^2).*(abs(psi_hat).^2) + abs(psi)/mu);
J = sum(0.5/N*((2*pi*(-N/2:N/2-1)/L)'.^2).*(abs(psi_hat).^2) + abs(psi)/mu);
% J = sum(abs(psi)/mu*sqrt(L/N));