function [psi, psi_hat, J] = cpw(lambda, r, mu, N, h, k, Psi_hat, max_iter)
% Algorithm 1: computation of the k-th CPW

% psi is the k-th BCPW, and psi_hat is the Fourier Transform of psi.

% Psi_hat = [psi1_hat psi2_hat ... psi_{k-1}_hat], 
% where psi1, psi2, ... psi_{k-1} are BCPWs of lower levels,
% psi1_hat, psi2_hat, ... psi_{k-1}_hat are their Fourier transforms.

% lambda and mu are defined the same as they are in the CPW paper.

% N is the (even) number of points on the domain [0, L].
% h is the number of points for the size of lattice w. 
% h divides N.
% [u, u_hat] = fourier_orthogonal(N,h,Psi_hat,k);
% [u, ~] = sopwbasis(N,h,k); u = ifftshift(u); u_hat = fftshift(fft(u));
u = zeros(N,1);
u(N/2,1) = 1; u_hat = fftshift(fft(u));
psi = zeros(N, 1);
b = zeros(N, 1);
P = cell(floor(N/(2*h))+1,1);
for omega=0:floor(N/(2*h))
    Psi_hat_omega = Psi_hat_matrix(Psi_hat, omega, N, h, k);
    n = size(Psi_hat_omega,1);
    P{omega+1,1} = eye(n)-Psi_hat_omega*Psi_hat_omega'/h;
end
% max_iter=2000; 
J = zeros(max_iter,1);
epsilon0=1e-6;
for jj=1:max_iter
    if jj==1 
        psi = u;
        psi_hat = u_hat;
    else
        [psi, psi_hat] = update_psi(u-b, N, h, k, lambda, r, P); 
    end
    psi = real(psi);       
    u = soft_shrink(psi+b, 1/(lambda*mu));
    b = b + psi - u;
    J(jj,1) = J_psi(mu, N, psi, psi_hat);
    if jj>1 && abs(J(jj,1)-J(jj-1,1))/J(jj-1,1)<epsilon0
        break;
    end
%     if jj>10 && J(jj,1)>J(jj-1,1)
%         break;
%     end
end
% psi = u;
% psi_hat = fft(u);
% norm(psi - u)/norm(u)