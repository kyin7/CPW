L = 16; w = 4;
N = 8192;
n = L/w;%total number of shifts
h = w/L*N; %should be an integer
mu = 5; lambda = 1e2*sqrt(L)/mu; r = 10*lambda;
alpha = 2; % mu_k = mu/k^alpha for level-k lambda and r should scale too (in function cpws)!
max_level = 1;

max_iter = 2000;

[psi, Psi_hat, J_Psi] = cpws(L, w, mu, alpha, lambda, r,  N, max_level, max_iter);
sgn_psi = zeros(N,1);
for l = 1:N
    if abs(psi(l)) < 0.000005
        sgn_psi(l) = 0;
    else
        sgn_psi(l) = sign(psi(l));
    end
end

dif_psi=(N/L)*diff(psi);%numerical differentiation, rescaling by mesh size 
dif2_psi=(N/L)*(N/L)*diff(psi,2);
PSI = zeros(n,N);% contains the shifts of psi
dif_PSI = zeros(n,N-1);% contains the derivatives of shifts of psi
lambda_vec =zeros(n,1);% lagrange multipliers
for l=1:n
    PSI(l,:) = circshift(psi,h*(l-1));
    dif_PSI(l,:) = circshift(dif_psi,h*(l-1));
    lambda_vec(l) = (L/N)*dot(dif_psi,dif_PSI(l,:)) + (1/mu)*(L/N)*dot(sgn_psi,PSI(l,:));
end

deviation = zeros(N-2,1);
for ind = 1:N-2
LHS = -dif2_psi(ind) + (1/mu)*sgn_psi(ind);
RHS = dot(lambda_vec, PSI(:,ind));
deviation(ind) = LHS - RHS;
end
figure;
plot((L/(N-2))*(0:N-3)-8, deviation);