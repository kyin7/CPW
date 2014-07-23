function psi2 = update_psi_V2(opt, level, rhs)
% solve the Schrodinger equation 
% -1/2*Laplace(psi) + (V+lambda+gamma).*psi = rhs
L = opt.L;
N = opt.N;
V = opt.V;
h = opt.h;
G = opt.G;
V_bar = opt.V_bar;
V_tilde = opt.V_tilde; % V_tilde = fft(V_bar);
lambda = opt.lambda; 
% lambda = opt.lambda/level^2; 
gamma = opt.gamma;
xx = opt.x -opt.center;
xx(xx<-L/2) = xx(xx<-L/2) + L;
xx(xx>L/2) = xx(xx>L/2) - L;
xx = abs(opt.x - opt.center);
tol = 1e-6;
maxit = 10;
% restart = 10; [psi, flag] = gmres(@Hfun, rhs, restart, tol, maxit, @Mfun);
% if flag~=0
%     disp('gmres not converged');
% end

% % solve the elliptic equation by Neumann Series Expansion
% rhs1 = Mfun(rhs);
% psi1_old = rhs1;
% for ii=1:maxit
%     psi1 = rhs1 - Mfun(V.*psi1_old);
%     if norm(psi1-psi1_old)/norm(psi1) < 1e-5
%         break;
%     else
%         psi1_old = psi1;
%     end
% end
% if ii==maxit
%     disp('Neumann Series Expansion not converged');
% end
% err1 = Hfun(psi1)-rhs;

% % solve by Neumann Series Expansion in k-space
% rhs_hat = fft(rhs);
% psi_hat = zeros(N,1);
% for k=0 : N/(2*h)
%     % solve sub-problems    
%     rhs_k = rhs_hat(k + G + 1); 
%     restart = 10; [psi_k, flag] = gmres(@H_bar, rhs_k, restart, tol, maxit, @M_bar);
%     if flag~=0
%         disp('gmres not converged');
%     end
%     % assemble psi_hat_k's into psi_hat
%     psi_hat(k + G + 1) = psi_k;
% end
% % fourier coefficients at negative wave numbers are conjugate of that of at
% % positive wave numbers
% for k= N/(2*h)+1 : N/h-1
%     psi_hat(k + G + 1) = conj(psi_hat(-k - G + N + 1));
% end
% % solution in the physical space
% psi2 = ifft((psi_hat)); psi2 = real(psi2);

% solve by Neumann Series Expansion in k-space
rhs_hat = fft(rhs);
psi_hat = zeros(N,1);
G2 = [0:h/2-1 -h/2:-1]'*N/h;
for k=0 : N/(2*h)
    % solve sub-problems
    k2 =  G2 + k;
    D_k = 0.5*(2*pi*k2/L).^2+lambda+gamma;
    rhs_k = rhs_hat(k + G + 1); 
    rhsk = rhs_k./D_k;
    psik_old = rhsk;   
    for ii=1:maxit
%         psi_k = rhsk - convolve(V_bar, psik_old)/N./D_k;
        psi_k = rhsk - convolution(psik_old)./D_k;
        if norm(psi_k-psik_old)/norm(psi_k) < tol
            break;
        else
            psik_old = psi_k;
        end
    end
    if ii==maxit
        disp('Neumann Series Expansion not converged');
    end
    % assemble psi_hat_k's into psi_hat
    psi_hat(k + G + 1) = psi_k;
end
% fourier coefficients at negative wave numbers are conjugate of that of at
% positive wave numbers
for k= N/(2*h)+1 : N/h-1
    psi_hat(k + G + 1) = conj(psi_hat(-k - G + N + 1));
end
% solution in the physical space
psi2 = ifft((psi_hat)); psi2 = real(psi2);


    function u = Hfun(f) % u = (-1/2*Laplace + V+lambda+gamma) f
        k2 = [0:N/2-1 -N/2:-1]';
        f1 = fft(f);
        f2 = real(ifft(((1i*2*pi/L*k2).^2.*f1))); % f2 = Laplace(f)
        u = -0.5*f2 + (V + 0e1*xx.^2 + lambda+gamma).*f;
    end

    function u = Mfun(f) % u = (-1/2*Laplace + lambda + gamma)^(-1) f
        k2 = [0:N/2-1 -N/2:-1]';
        f1 = fft(f);
        f2 = f1./(-0.5*(1i*2*pi/L*k2).^2+lambda+gamma);
        u = real(ifft(f2));
    end

    function u = H_bar(f)
        k2 = [0:h/2-1 -h/2:-1]'*N/h + k;
        D_k = 0.5*(2*pi*k2/L).^2+lambda+gamma;
        u = D_k.*f + convolve(V_bar, f)/N;
    end
    
    function u = M_bar(f)
        k2 = [0:h/2-1 -h/2:-1]'*N/h + k;
        u = f./(0.5*(2*pi*k2/L).^2+lambda+gamma);
    end

    function u = convolution(f)
        u = ifft(V_tilde.*fft(f))/N;
    end
end
