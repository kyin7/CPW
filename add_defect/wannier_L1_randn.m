function [psi, psi_hat, J, J_min] = wannier_L1_randn(k, opt, Psi_hat)
% Algorithm 1: computation of the k-th CPW
%
% psi is the k-th BCPW, and psi_hat is the Fourier Transform of psi.
%
% Psi_hat = [psi1_hat psi2_hat ... psi_{k-1}_hat],
% where psi1, psi2, ... psi_{k-1} are BCPWs of lower levels,
% psi1_hat, psi2_hat, ... psi_{k-1}_hat are their Fourier transforms.
%
% lambda is the parameter for shift-orthogonal constraint.
% N(even) is the number of points on the domain [0, L].
% h is the number of points for the size of lattice w.
% h divides N.
% use a random vector as the initial guess
N = opt.N;
n = opt.n;
h = opt.h;
% L = opt.L;
max_iter = opt.max_iter;
lambda = opt.lambda;
gamma = opt.gamma;
%

% projection matrix for shift-orthogonality
P = cell(floor(N/(2*h))+1,1);
for omega=0:floor(N/(2*h))
    P{omega+1,1} = Psi_hat_matrix(Psi_hat, omega, N, h, k);
    %     Psi_hat_omega = Psi_hat_matrix(Psi_hat, omega, N, h, k);
    %     n = size(Psi_hat_omega,1);
    %     P{omega+1,1} = eye(n)-Psi_hat_omega*Psi_hat_omega'/h;
end

J = zeros(max_iter,1);
J_min = 0;
epsilon0=1e-4;

np = 2*n; % try np different initial data in order to get global min
for ii=1:np
    u = zeros(N,1); u(N/2+(ii-n)*N/n/2,1) = 1;
%     u = randn(N,1);
%     u = exp(-(opt.x-opt.L/2).^2) + 1e-2*randn(N,1); 
%     u = u/norm(u);
%     opt.center = find_center(u, opt.x);
    v = u;
    b = zeros(N, 1);
    c = zeros(N, 1);
    psi1 = zeros(N, 1);
    psi2 = u;
    J2 = zeros(max_iter,1);
    for jj=1:max_iter
        % first minimize the augmented lagrangian wrt psi
        % psi = argmin <psi, H psi> + lambda/2*|psi-u+b|^2 +
        % gamma/2*|psi-v+c|^2
%         rhs = lambda*(u-b);
        rhs = lambda*(u-b) + gamma*(v-c);
        % solve (H+lambda+gamma) psi = rhs
%         psi2 = update_psi_V2(opt, k, rhs);
%         psi2 = update_psi_V(opt, k, rhs);
        psi2 = update_psi(opt, k, rhs);
        % secondly project psi+b onto shift-orthogonal set
        [u, u_hat] = proj_sor(psi2+b, k, opt, P);
        % L1 minimization by soft-thresholding 1/(gamma*mu)
        v = L1_step(psi2+c, opt);
        % make sure psi is a real function
        psi2 = real(psi2);
%         opt.center = find_center(psi2, opt.x);
%         opt.center = dot(u, opt.x.*u);
        % update auxillary variable(s)
        b = b + psi2 - u;
        c = c + psi2 - v;
        % evaluate energy
        J2(jj,1) = J_psi(opt, u, u_hat);
        % stopping criteria
        if jj>1 && norm(psi1-psi2)<epsilon0 %&& abs(J2(jj,1)-J2(jj-1,1))/abs(J2(jj-1,1))<epsilon0
            break;
        end
        psi1 = psi2;
%         if jj>1 && (norm(psi2-u) < epsilon0) && abs(J2(jj,1)-J2(jj-1,1))/abs(J2(jj-1,1))<epsilon0
%             break;
%         end
    end
    if ii==1 || J2(jj,1) < J_min
        psi = psi2;
        J_min = J2(jj,1);
        J = J2;
%         str = ['Level ' num2str(k) ': current energy = ' num2str(J_min)];
%         disp(str);
    end
end
J_min = J_min - sum(abs(psi)/opt.mu) * sqrt(opt.L/opt.N);
str = ['Level ' num2str(k) ': minimum energy = ' num2str(J_min)];
disp(str);
% if k<opt.max_level
%     center_point = find_center(psi, opt.x);
%     psi = shift(psi, floor(center_point)-opt.L/2, h);
% end
psi_hat = fftshift(fft(psi));
% norm(psi - u)/norm(u)
