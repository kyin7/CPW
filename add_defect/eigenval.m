function eigenvalues = eigenval(opt, max_number)
% eigenvalues of Hamiltonian, in acsending order
N = opt.N;
V = opt.V;
L = opt.L;
opt.isreal = true;
max_eigval = abs(eigs(@afun, N, 1, 'LM', opt));
eigenvalues = eigs(@bfun, N, max_number, 'LM', opt) + max_eigval;
    function u = afun(f) % u = -1/2*Laplace(f) + (V+lambda).*f
        G = [0:N/2-1 -N/2:-1]';
        f1 = fft(f);
        f2 = real(ifft(((1i*2*pi/L*G).^2.*f1))); % f2 = Laplace(f)
        u = -0.5*f2 + V.*f;
    end
    function u = bfun(f) % u = -1/2*Laplace(f) + (V+lambda).*f
        G = [0:N/2-1 -N/2:-1]';
        f1 = fft(f);
        f2 = real(ifft(((1i*2*pi/L*G).^2.*f1))); % f2 = Laplace(f)
        u = -0.5*f2 + (V-max_eigval).*f;
    end
figure; plot(sort(eigenvalues), '*');
% V = opt.V(1:n,1);
% L = opt.L;
% lambda = opt.lambda;
% gamma = opt.gamma;
% max_number = n;
% eigenvalues = eigs(@bfun, n, max_number)-lambda-gamma;
%     function u = afun(f) % u = -1/2*Laplace(f) + (V+lambda).*f
%         G = [0:n/2-1 -n/2:-1]';
%         f1 = fft(f);
%         f2 = real(ifft(((1i*2*pi/L*G).^2.*f1))); % f2 = Laplace(f)
%         u = -0.5*f2 + (V+lambda+gamma).*f;
%     end
%     function v = bfun(f_hat)
%         G = [0:n/2-1 -n/2:-1]'+k;
%         f = ifft(f_hat);
%         v = -0.5*(1i*2*pi/L*G).^2.*f_hat + fft((V+lambda+gamma).*f);
%     end
% % figure;plot(real(eigenvalues(end:-1:end-79)), '*');
% figure;plot(real(eigenvalues(end:-1:1)));
end

% %%
% N = opt.N;
% V = opt.V*0;
% L = opt.L;
% x = [0:N/2-1 -N/2:-1]'*L/N;
% f = sin(2*pi*x/L);
% k = [0:N/2-1 -N/2:-1]';
% f1 = fft(f);
% f2 = real(ifft(((1i*2*pi/L*k).^2.*f1))); % f2 = Laplace(f)
% figure;plot(f2);