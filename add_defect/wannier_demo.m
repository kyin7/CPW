function [Psi, Psi_hat, J_Psi, Psi_energy, opt]= wannier_demo
%% parameters
n = 12; % number of unit cells in a supercell;
w = 1; L = n*w; 
h = 200;
N = h*n;
G = (0:h-1)'*n; % reciprocal lattice

% d = 1; s=1/128;
% L = L*s; w = w*s;

V0 = 1e2; % maximum depth of the potential wells

% No L1
% mu = 1e10; gamma = 0;% no L1
% lambda = 10*V0; % coefficient for shift-orthogonality


% add L1
% mu = 10/(sqrt(L));
% gamma = 1e3; % coefficient for L1 minimization
% lambda = 1e3; % coefficient for shift-orthogonality

mu = 2e1/(V0);
gamma = 0.5e0*V0; % coefficient for L1 minimization
lambda = 5*V0; % coefficient for shift-orthogonality



max_iter = 1000;
max_level = 4;

opt.L = L;
opt.N = N;
opt.w = w;
opt.h = h; %% the size of the shift is multiple of unit cell
opt.G = G;
opt.mu = mu;
opt.lambda = lambda;
opt.gamma = gamma;
opt.max_iter = max_iter;
opt.max_level = max_level;
% opt.center = L/2;
x = ((0:N-1)*L/N)';
opt.x = x;
opt.n = n;
%% potential
V = zeros(N,1);
np = 2*L/w-1;
xp = (0:np+1)*L/(np+1);
% V0 = 1e3; 
sigma = w/10;
for ii=1:np+2
    if mod(ii,2)==1
        rr = 1;
    else 
        rr = 0.6;
    end
    V = V - rr*V0*exp(-(x-xp(ii)).^2/(2*sigma^2));
end
V(1:h) = V(h+1:2*h);
V(end-h:end) = V(end-2*h:end-h);
% shift the potential V;
V = circshift(V, h/2 * 1);
V_hat = fft(V);
V_bar = V_hat(G+1);
V_tilde = fft(V_bar);
opt.V = V;
opt.V_bar = V_bar;
opt.V_tilde = V_tilde;

% add defect to the potential
V_d = zeros(N,1);
for ii=1:np+2
    if ii==floor((np+1)/2) % add the defect here
        rr = 1.5;
    elseif mod(ii,2)==1
        rr = 1;
    else 
        rr = 0.6;
    end
    V_d = V_d - rr*V0*exp(-(x-xp(ii)).^2/(2*sigma^2));
end
% V_d(1:h/2) = V_d(h+1:h+h/2);
% V_d(end-h/2:end) = V_d(end-h-h/2:end-h);
% shift the potential V;
V_d = circshift(V_d, h/2*1);
V_d_hat = fft(V_d);
V_d_bar = V_d_hat(G+1);
V_d_tilde = fft(V_d_bar);
opt.V_d = V_d;
opt.V_d_bar = V_d_bar;
opt.V_d_tilde = V_d_tilde;

% opt.V = V_d;
% opt.V_bar = V_d_bar;
% opt.V_tilde = V_d_tilde;

% opt.lambda = 10*sqrt(h)*norm(opt.V_bar);
% close all
% figure;
% plot((0:N-1)*L/N,V, 'LineWidth', 2);
% xlim([0,L]);
% % title(['Kronig-Penny potential: ' 'L = ' num2str(L) ' w = ' num2str(w) ],'FontSize',15)
% set(gca, 'XTick', 0:w:L, 'FontSize', 15);
% grid on;
%% Wannier Functions
% [Psi, Psi_hat, J_Psi] = cpws(opt, max_level);
[Psi, Psi_hat, J_Psi, Psi_energy] = wannier(opt, max_level);
% [Psi, Psi_hat, J_Psi] = wannier_sopw(opt, max_level);

% plot the result
% plot_first_8(Psi, Psi_energy, opt);
% plot_first_4(Psi, Psi_energy, opt);
% plot potential
% figure; plot((0:N-1)*L/N,V_d, 'LineWidth', 2);
% xlim([0,L]);
% ylim([min(V_d), 0]);
% set(gca, 'XTick', 0:w:L, 'FontSize', 15);

% plotresult(Psi, Psi_hat, J_Psi(20:end,:), opt);

%% scaling formula
% d = 1; s=0.4;
% L = L*s; w = w*s;
% mu = mu*s^(2+d/2);
% 
% max_level = 4;
% tic
% [Psi_s, Psi_s_hat, J_Psi_s] = cpws(L, w, mu, lambda, r,  N, max_level, max_iter);
% toc
% Psi_s = Psi_s * s^(d/2);
% Psi_s_hat = Psi_s_hat * s^(d/2);
% J_Psi_s = J_Psi_s * s^(3/2+d/2);
% plotresult(Psi_s, Psi_s_hat, J_Psi_s, L, w, mu);

%%
% Psi_hat_sq = zeros(N,1);
% for ii=1:N
%     Psi_hat_sq(ii,1) = norm(Psi_hat(ii,:))^2;
% end
% figure;plot(-N/2:N/2-1, Psi_hat_sq);
% ylim([0, sqrt(N)+1]);
% title('sum of absolute value squared of Fourier coefficients of all levels', 'FontSize', 12);