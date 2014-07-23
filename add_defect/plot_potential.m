function plot_potential(opt)
N = opt.N;
V = opt.V;
L = opt.L;
w = opt.w;
figure;
plot((0:N-1)*L/N,V, 'LineWidth', 2);
xlim([0,L]);
ylim([min(V), 0]);
% title(['Kronig-Penny potential: ' 'L = ' num2str(L) ' w = ' num2str(w) ],'FontSize',15)
set(gca, 'XTick', 0:w:L, 'FontSize', 15);
grid on;