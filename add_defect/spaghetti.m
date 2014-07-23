function [eigens_k, eigs_k] = spaghetti(opt, Psi, Psi_hat)
% eigens_k are eigenvalues calculated by Fourier basis
% eigs_k are eigenvalues calculated by Wannier modes
l = 5;
m = size(Psi_hat, 2);
levels = 1:m;
eigens_k = zeros(l*opt.N/opt.h/2+1, m);
eigs_k = zeros(l*opt.N/opt.h/2+1, m);
for kk=0:l*opt.N/opt.h/2
    % dispersion in plane wave basis
    eigens_kk = k_dispersion(opt, kk, l, m);
    eigens_k(kk+1,levels) = eigens_kk(levels);
    % dispersion in wannier basis
    eigs_kk = k_wannier(Psi, Psi_hat, opt, kk, l);
    eigs_k(kk+1,levels) = eigs_kk(levels);
end

figure;
for level = levels
    plot(eigens_k(:,level),'-r','LineWidth', 2); hold on;
%     hold all;
end
for level = levels
end
% 
% figure
for level = levels
    plot(eigs_k(:,level),'*b','LineWidth', 1); hold on;
%     hold all;
end
xlim([1,l*opt.N/opt.h/2+1]);
xlabel('k', 'FontSize', 15)
ylabel('E(k)', 'FontSize', 15);
set(gca,'XTick',[1 l*opt.N/opt.h/2+1]);
% set(gca, 'XTickLabel', {'0', 'p/a'}, 'FontSize', 15, 'fontname','symbol');
set(gca, 'XTickLabel', []); 
% text(1, -115, '0', 'FontSize', 15);
% text(l*opt.N/opt.h/2+1, -115, '\pi/a', 'FontSize', 15);
set(gca, 'FontSize', 15);

