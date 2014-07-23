function [u, u_hat] = proj_sor(v, k, opt, P)
% u = min||u-v||^2 s.t. <u(.), u(.+lh)>=delta(l) and <u, psi_j>=0 for j<k
% N is an even number, and h divides N.
% psi_j is the j-th CPW
N = opt.N;
h = opt.h;
L = opt.L;
v_hat = fftshift(fft(v));
u_hat = zeros(N,1);
p_shift = N/2+1;

% omega is zero
omega = 0;
x_p = proj_sor_sub(v_hat, k, opt, omega, P{omega+1,1});
if mod(h,2)==0
    p = (-floor(h/2) : 0)*N/h;
%     if norm(x_p)<1e-6*sqrt(h)
%         x_p(floor(h/2)+1,1)=sqrt(h);
%     end
else
    p = (-floor(h/2-1/2) : 0)*N/h;
%     if norm(x_p)<1e-6*sqrt(h/2)
%         x_p(floor(h/2-1/2)+1,1)=sqrt(h);
%     end
end

u_hat(p + p_shift,1) = x_p;
p_prime = (1:floor(h/2-1/2))*N/h;
u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift,1));

% omega is non-zero
for omega = 1:floor(N/(2*h)-1/2)
    p = omega+(-floor(h/2):floor(h/2-1/2))*N/h;
    x_p = proj_sor_sub(v_hat, k, opt, omega, P{omega+1,1});
%     if norm(x_p)<1e-6*sqrt(h)
%         x_p(floor(h/2)+1,1)=sqrt(h);
%     end
    u_hat(p + p_shift,1) = x_p;
    u_hat(-p + p_shift,1) = conj(u_hat(p + p_shift,1));
end

if floor(N/(2*h)) == N/(2*h)
    omega = N/(2*h);
    x_p = proj_sor_sub(v_hat, k, opt, omega, P{omega+1,1});
    if mod(h,2)==1
        p = -N/2 : N/h : -N/(2*h);
    else
        p = -N/2 + N/(2*h) : N/h : -N/(2*h);   
    end
%     if norm(x_p)<1e-6*sqrt(h)
%         x_p(floor(h/2),1)=sqrt(h/2);
%     end
    u_hat(p + p_shift,1) = x_p;
    p_prime = N/(2*h) + (0:floor(h/2)-1)*N/h;
    u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift));
end    
% u_hat = u_hat * sqrt(N/L);
% inverse fft
u = ifft(ifftshift(u_hat));
u = u/norm(u);
u_hat = fftshift(fft(u));