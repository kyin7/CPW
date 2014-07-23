function psi = update_psi_V(opt, level, rhs)
% solve the Schrodinger equation 
% -1/2*Laplace(psi) + (V+lambda+gamma).*psi = rhs
L = opt.L;
N = opt.N;
% if level<opt.max_level
    V = opt.V;
% else
%     V = opt.V_d;
% end
lambda = opt.lambda; 
% lambda = opt.lambda/level^2; 
gamma = opt.gamma;
% rhs = lambda*(u-b);
xx = opt.x -opt.center;
xx(xx<-L/2) = xx(xx<-L/2) + L;
xx(xx>L/2) = xx(xx>L/2) - L;
% xx = abs(opt.x - opt.center);
tol = 1e-4;
maxit = 10;
% k = [0:N/2-1 -N/2:-1]';
% rhs_hat = fft(rhs);
% rhs_f = 0.5 * real(ifft(((2*pi/L*k).^2.*rhs_hat))); % rhs_f = - 1/2 * Laplace(rhs)
% D = dot(rhs, real(rhs_f))/dot(rhs, rhs);
% [psi, flag] = pcg(@afun, rhs, tol, maxit, @mfun);
restart = 10; [psi, flag] = gmres(@afun, rhs, restart, tol, maxit, @mfun);
if flag~=0
    disp('gmres not converged');
end

    function u = afun(f) % u = -1/2*Laplace(f) + (V+lambda).*f
        k = [0:N/2-1 -N/2:-1]';
        f1 = fft(f);
        f2 = real(ifft(((1i*2*pi/L*k).^2.*f1))); % f2 = Laplace(f)
        u = -0.5*f2 + (V+0e2*xx.^2+lambda+gamma).*f;
    end

    function u = mfun(f)
        % if the potential energy is dominant over the kinetic energy,
        % we use the kinetic part of H as the preconditioner     
%         u = f./(D + V+lambda+gamma);
%         u = f./(V+lambda+gamma);
%         u = f./V;
        % if the kinetic energy is dominant, we use the solution operator
        % of Poisson equation as the preconditioner
        k = [0:N/2-1 -N/2:-1]';
        f1 = fft(f);
        f2 = f1./(-0.5*(1i*2*pi/L*k).^2+lambda+gamma);
        u = real(ifft(f2));
    end
end
