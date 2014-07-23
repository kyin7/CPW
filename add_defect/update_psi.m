function psi1 = update_psi(opt, level, rhs)
% solve the Schrodinger equation 
% -1/2*Laplace(psi) + (V+lambda+gamma).*psi = rhs
N = opt.N;
L = opt.L;
% if level<opt.max_level
    V = opt.V;
% else
%     V = opt.V_d;
% end
lambda = opt.lambda; 
gamma = opt.gamma;
tol = 1e-5;
maxit = 10;

% solve the elliptic equation by Neumann Series Expansion
rhs1 = Mfun(rhs);
psi1_old = rhs1;
for ii=1:maxit
    psi1 = rhs1 - Mfun(V.*psi1_old);
    if norm(psi1-psi1_old)/norm(psi1) < tol
        break;
    else
        psi1_old = psi1;
    end
end
if ii==maxit
    disp('Neumann Series Expansion not converged');
end

    function u = Mfun(f) % u = (-1/2*Laplace + lambda + gamma)^(-1) f
        f1 = fft(f);
        k2 = [0:N/2-1 -N/2:-1]';
        D = (-0.5*(1i*2*pi/L*k2).^2+lambda+gamma);
        f2 = f1./D;
        u = real(ifft(f2));
    end
%     function u = Mfun(f) % u = (-1/2*Laplace + lambda + gamma)^(-1) f
%         f1 = fft(f);
%         dx = L/N;
%         C = zeros(N,1);
%         C(1,1) = 1/dx^2 + lambda + gamma;
%         C(2,1) = -0.5/dx^2;
%         C(N,1) = -0.5/dx^2;
%         D = fft(C);
%         f2 = f1./D;
%         u = real(ifft(f2));
%     end
end