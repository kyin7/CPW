function plot_first_8(Psi, Psi_energy, opt)
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

% figure;
% set(gcf,'Position', [379   405   433   555]);
% set(gcf,'PaperPositionMode','auto');
for jj=5:8
    subplot(4,2,jj);
    psi = Psi(:,jj);
    plot(L/N*(0:N-1), psi, 'LineWidth', 2);
    set(gca, 'XTick', 0:w:L, 'FontSize',15);
    grid on;
%     xlim([0,L]);
    xlim([0,L]);
    ylim([-2.5,2.5]);
%     xlabel('x');
%     ylabel('\psi(x)');
    title(['  Level-' num2str(jj) ' Energy = ' num2str(Psi_energy(jj))],'FontSize',15);
end

figure;
for jj=9:12
    subplot(4,2,jj-8);
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

for jj=13:16
    subplot(4,2,jj-8);
    psi = Psi(:,jj);
    plot(L/N*(0:N-1), psi, 'LineWidth', 2);
    set(gca, 'XTick', 0:w:L, 'FontSize',15);
    grid on;
%     xlim([0,L]);
    xlim([0,L]);
    ylim([-2.5,2.5]);
%     xlabel('x');
%     ylabel('\psi(x)');
    title(['  Level-' num2str(jj) ' Energy = ' num2str(Psi_energy(jj))],'FontSize',15);
end

figure;
for jj=17:20
    subplot(4,2,jj-16);
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

for jj=21:24
    subplot(4,2,jj-16);
    psi = Psi(:,jj);
    plot(L/N*(0:N-1), psi, 'LineWidth', 2);
    set(gca, 'XTick', 0:w:L, 'FontSize',15);
    grid on;
%     xlim([0,L]);
    xlim([0,L]);
    ylim([-2.5,2.5]);
%     xlabel('x');
%     ylabel('\psi(x)');
    title(['  Level-' num2str(jj) ' Energy = ' num2str(Psi_energy(jj))],'FontSize',15);
end