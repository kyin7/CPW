function plot_result(psi, psi_hat, J, L, w, mu, alpha, N, k)
% plot_result plot the BCPW psi, absolute value squared of the Fourier
% Transform |psi_hat|^2, total energy evolution J_psi
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/2])
% psi
subplot(2,2,1);
plot(L/N*(0:N-1), psi);
xlabel('x');
ylabel('\psi(x)');
title(['L = ' num2str(L) ' w = ' num2str(w) ' \mu = ' num2str(mu/k^alpha) '  BCPW-' num2str(k)],'FontSize',15);
% |psi_hat|^2
subplot(2,2,2);
plot(-N/2:N/2-1,abs(psi_hat).^2);
xlabel('\omega');
text('Interpreter','latex',  'String', '$|\hat{\psi}(\omega)|^2$', 'Position',[2*L/(w) N^2/L^2*w/3]);
title('Fourier Analysis', 'FontSize', 15);
% J(psi)
subplot(2,2,3);
plot(find(J>0), (2*pi/L)^2 * J(J>0));
total_ind = find(J>0);
max_ind = total_ind(end);
title(['Total Energy: < \psi, H \psi > + 1/mu |\psi |_1 = ' num2str((2*pi/L)^2 * J(max_ind)) ],'FontSize',15);
xlabel('iterations');
ylabel('J(\psi)');