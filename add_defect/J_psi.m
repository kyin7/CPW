function J = J_psi(opt, psi, psi_hat)
% J(psi) = <psi, H psi> + 1/mu |psi|_1
V = opt.V;
N = opt.N;
L = opt.L;
% J = 0.5*(N/L)^2*norm(diff(psi))^2 + dot(psi, V.*psi);
k = (-N/2:N/2-1)';
J = -0.5*sum((1i*2*pi*k/L).^2.*(abs(psi_hat).^2))/N + dot(psi, V.*psi);
% J = - N/L*0.5/L*sum((1i*2*pi*k/N).^2.*(abs(psi_hat).^2)) + dot(psi, V.*psi);
if opt.mu > 0
    J = J + sum(abs(psi)/opt.mu) * sqrt(L/N);
end

% k = [0:N/2-1 -N/2:-1]';
% a1 = dot(V, ifft((1i*2*pi*k/N).^2.*fft(V)));
% a2 = 1/N*dot(fft(V), (1i*2*pi*k/N).^2.*fft(V));