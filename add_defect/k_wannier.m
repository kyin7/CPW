function [eigs_k] = k_wannier(Psi, Psi_hat, opt, k, l)
% h_{ii,jj}(k) = \sum_p <psi^{ii}(x+p*w), H^ psi^{jj}(x)> e^{-i k p*w}
% eigs_k are eigenvalues of h

m = size(Psi_hat,2);
% Psi_hat_k = Psi_hat(omega, :);
N = opt.N;
h = opt.h;
V = repmat(opt.V, l ,1);
L = l*opt.L;
n = l*N;

h_k = zeros(m);
psi_ii = zeros(n,1);
psi_jj = zeros(n,1);
omega = k + (0:h-1)*n/h + 1;
for ii=1:m
    for jj=1:m
%         psi = Psi(:,jj);
%         h_psi_hat = fftshift(fft(Hfun(psi)));
        psi_ii(1:N,1) = Psi(:,ii);
        psi_jj(1:N,1) = Psi(:,jj);
        psi_ii_hat = fft(psi_ii);
        psi_jj_hat = fft(psi_jj);
        h_psi_hat = hfun(psi_jj, psi_jj_hat);
        h_k(ii,jj) = L/(n*h) * dot(psi_ii_hat(omega, 1), h_psi_hat(omega, 1));
    end
end
eigs_k = sort(real(eig(h_k)));

% H_k = zeros(m);
% psi_ii = zeros(n,1);
% psi_jj = zeros(n,1);
% for ii=1:m
%     for jj=1:m
%         psi_ii(1:N,1) = Psi(:,ii);
%         psi_jj(1:N,1) = Psi(:,jj);
%         psi_ij = zeros(n/h,1);
%         H_psi_jj = Hfun(psi_jj);
%         for kk=0:n/h-1
%             psi_ij(kk+1,1) = L/n*dot(shift(psi_ii, kk, h), H_psi_jj);
%         end
%         psi_ij_hat = fft(psi_ij);
%         H_k(ii,jj) = psi_ij_hat(k+1);
%     end
% end
% eigs_k = sort(real(eig(H_k)));
%     

% figure;plot(eigs_k, '*');

% subroutines
    function v = hfun(f, f_hat)
        G = [0:n/2-1 -n/2:-1]';
%         f = ifft(f_hat);
        v = (-0.5*(1i*2*pi/L*G).^2.*f_hat) + fft(V.*f);
    end

%     function u = Hfun(f)
%         G = [0:n/2-1 -n/2:-1]';
%         f_hat = fft(f);
%         u = ifft(-0.5*(1i*2*pi/L*G).^2.*f_hat) + V.*f;
%     end
% end subroutines
end
