function plot_first_4(Psi, L, w, mu, alpha)
N = size(Psi,1);
figure;
for jj=1:4
    subplot(2,2,jj);
    psi = Psi(:,jj);
    plot(L/N*(0:N-1), psi);
    xlabel('x');
    ylabel('\psi(x)');
    axis('tight');
    title(['L = ' num2str(L) ' w = ' num2str(w) ' \mu = ' num2str(mu/jj^alpha) '  BCPW-' num2str(jj)],'FontSize',15);
end
figure;
for jj=1:4
    subplot(2,2,jj);
    psi = Psi(:,jj+4);
    plot(L/N*(0:N-1), psi);
    xlabel('x');
    ylabel('\psi(x)');
    axis('tight');
    title(['L = ' num2str(L) ' w = ' num2str(w) ' \mu = ' num2str(mu/(jj+4)^alpha) '  BCPW-' num2str(jj+4)],'FontSize',15);
end

