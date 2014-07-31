function [Psi, Psi_hat, J_Psi, L, w, N, h, mu, alpha]= cpw_demo
clear all;
%% parameters
% L = 64*2; w = 8;
% N = 256*4;  
% h = w/L*N; %should be an integer
% mu = 200;
% d = 1; s=1/8;
% L = L*s; w = w*s;
% mu = mu*s^(2+d/2);

w = 1;
L = 6*w;
h = 100;
N = h*L/w;
mu = 0.5;

lambda = 1e2*sqrt(L)/mu; r = 10*lambda;
alpha = 0; % mu_k = mu/k^alpha for level-k lambda and r should scale too (in function cpws)!
max_level = 1;

max_iter = 2000;

% opt.L = L;
% opt.N = N;
% opt.alpha = alpha;
% opt.h = h;
% opt.mu = mu;

%% BCPWs
[Psi, Psi_hat, J_Psi] = cpws(L, w, mu, alpha, lambda, r,  N, max_level, max_iter);
% plot the result
% close all
% plot_first_4(Psi, L, w, mu, alpha);
plotresult(Psi, Psi_hat, J_Psi, L, w, mu, alpha);

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