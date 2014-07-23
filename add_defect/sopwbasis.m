function [u, u_hat] = sopwbasis(N, h, k)
% u is the k-th shift-orthogonal basis, whose length is N, and
%   <u(x-lh), u(x)> = delta(l)
% N is even and divisible by h

u_hat = zeros(N,1);
omega_offset = N/2 + 1;
if floor(N/(2*h)) == N/(2*h)
    if k==1
        omega = (k-1)*N/(2*h) : k*N/(2*h);
        u_hat(omega(1) + omega_offset) = sqrt(h);
        u_hat(omega(N/(2*h)+1) + omega_offset) = sqrt(h/2);
    elseif k==h
        omega = (k-1)*N/(2*h) : k*N/(2*h)-1;
        u_hat(omega(1) + omega_offset) = 1i^(k-1)*sqrt(h/2);
        u_hat(1) = sqrt(h);        
    else
        omega = (k-1)*N/(2*h) : k*N/(2*h);
        u_hat(omega(1) + omega_offset) = 1i^(k-1)*sqrt(h/2);
        u_hat(omega(N/(2*h)+1) + omega_offset) = 1i^(k-1)*sqrt(h/2);        
    end
    u_hat(omega(2:N/(2*h)) + omega_offset) = 1i^(k-1)*sqrt(h);
    u_hat(-omega + omega_offset) = conj(u_hat(omega + omega_offset));
else % floor(N/(2*h)) < N/(2*h) 
    if k==1
        omega = 0 : (N-h)/(2*h);
        u_hat(omega + omega_offset) = sqrt(h);     
    elseif k==h
        omega = ((h-1)*N+h)/(2*h) : N/2-1;
        u_hat(omega + omega_offset) = 1i^(k-1)*sqrt(h);
        u_hat(1) = sqrt(h);
    elseif mod(k,2)==1 % k is odd
        omega = (k-1)*N/(2*h) : (k*N-h)/(2*h);
        u_hat(omega(1) + omega_offset) = 1i^(k-1)*sqrt(h/2);
        u_hat(omega(2:end) + omega_offset) = 1i^(k-1)*sqrt(h);
    else
        omega = ((k-1)*N+h)/(2*h) : k*N/(2*h);
        u_hat(omega(end) + omega_offset) = 1i^(k-1)*sqrt(h/2);
        u_hat(omega(1:end-1) + omega_offset) = 1i^(k-1)*sqrt(h);
    end 
    u_hat(-omega + omega_offset) = conj(u_hat(omega + omega_offset));
end
u = (ifft(fftshift(u_hat)));
