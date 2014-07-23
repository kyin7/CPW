function plot_first_4(Psi, Psi_energy, opt)
% N = opt.N;
% L = opt.L;
% w = opt.w;
% figure;
% set(gcf,'Position', [379   405   433   555]);
% set(gcf,'PaperPositionMode','auto');
% for jj=1:4
%     subplot(4,1,jj);
%     psi = Psi(:,jj);
%     plot(L/N*(0:N-1), psi(1:N), 'LineWidth', 2);
%     set(gca, 'XTick', 0:w:L, 'FontSize',15);
%     grid on;
%     xlim([0,L]);
% %     xlabel('x');
% %     ylabel('\psi(x)');
%     title(['  Level-' num2str(jj)],'FontSize',15);
% end

% set(gcf,'Position', [379   405   433   555]);
% set(gcf,'PaperPositionMode','auto');

% figure;
% for jj=5:8
%     subplot(4,1,jj-4);
%     psi = Psi(:,jj);
%     plot(L/N*(0:N-1), psi, 'LineWidth', 2);
%     set(gca, 'XTick', 0:w:L, 'FontSize',15);
%     grid on;
%     xlim([0,L]);
% %     xlabel('x');
% %     ylabel('\psi(x)');
%     title(['  Level-' num2str(jj)],'FontSize',15);
% end

N = opt.N;
L = opt.L;
w = opt.w;
figure;
% set(gcf,'Position', [379   405   433   555]);
% set(gcf,'PaperPositionMode','auto');
for jj=1:4
    subplot(4,2,jj);
    psi = Psi(:,jj);
    plot(L/N*(0:N-1), psi(1:N), 'LineWidth', 2);
    set(gca, 'XTick', 0:w:L, 'FontSize',15);
    grid on;
%     xlim([0,L]);
    xlim([0,L]);
    ylim([-2.5,2.5]);
%     xlabel('x');
%     ylabel('\psi(x)');
    title(['  Level-' num2str(jj) ' Energy = ' num2str(Psi_energy(jj))],'FontSize',15);
end