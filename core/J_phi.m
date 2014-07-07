function J = J_phi(L, mu, N, phi)
% J(phi) = <phi, H phi> + 1/mu*|phi|_1
phi = phi/norm(phi);
phi_hat = fftshift(fft(phi));
scaling = (2*pi/L)^2*sqrt(N/L);
mu_s = mu * scaling;
J = (2*pi/L)^2 * sum(0.5/N*((-N/2:N/2-1)'.^2).*(abs(phi_hat).^2) + abs(phi)/mu_s);