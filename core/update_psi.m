function [psi, psi_hat] = update_psi(v, N, h, k, lambda, r, P)
% psi = argmin <psi, H psi> + lambda/2*||psi-v||^2 s.t. <psi(.),
% psi(.+lh)>=delta(l) and <psi, psi_i>=0 for i<k
% N is an even number, and h divides N.
% psi is the k-th CPW
v_hat = fftshift(fft(v));
psi_hat = zeros(N,1);
p_shift = N/2+1;

% omega is zero
omega = 0;
x_p = solve_x(v_hat, lambda, r, N, h, k, omega, P{omega+1,1});
if mod(h,2)==0
    p = (-floor(h/2):0)*N/h;
else
    p = (-floor(h/2-1/2) : 0)*N/h;
end
psi_hat(p + p_shift,1) = x_p;
p_prime = (1:floor(h/2-1/2))*N/h;
psi_hat(p_prime + p_shift,1) = conj(psi_hat(-p_prime + p_shift,1));

% omega is non-zero
for omega = 1:floor(N/(2*h)-1/2)
    p = omega+(-floor(h/2):floor(h/2-1/2))*N/h;
    x_p = solve_x(v_hat, lambda, r, N, h, k, omega, P{omega+1,1});
    psi_hat(p + p_shift,1) = x_p;
    psi_hat(-p + p_shift,1) = conj(psi_hat(p + p_shift,1));
end

if floor(N/(2*h)) == N/(2*h)
    omega = N/(2*h);
    x_p = solve_x(v_hat, lambda, r, N, h, k, omega, P{omega+1,1});
    if mod(h,2)==1
        p = -N/2 : N/h : -N/(2*h);
    else
        p = -N/2 + N/(2*h) : N/h : -N/(2*h);   
    end
    psi_hat(p + p_shift,1) = x_p;
    p_prime = N/(2*h) + (0:floor(h/2)-1)*N/h;
    psi_hat(p_prime + p_shift,1) = conj(psi_hat(-p_prime + p_shift));
end    

% inverse fft
psi = ifft(ifftshift(psi_hat));
% psi = real(psi);
% ind_imag = abs(imag(psi_hat))<1e-6;
% psi_hat(ind_imag) = real(psi_hat(ind_imag));