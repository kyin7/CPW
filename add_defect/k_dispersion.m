function eigens = k_dispersion(opt, k, l, max_number)
% eigenvalues of k-dispersion matrix 
%
% H(k+G, k+G') = 1/2*(k+G)^2*delta_{G,G'} + V(G-G')
% 0<=k<opt.h/2
% V(.) is the Fourier transform of the potential
N = opt.N;
h = opt.h;
n = l*N;
if k > n/h
    error('k should be less than n/h');
elseif max_number > h
    error('max_number should be less than h');
end
V = opt.V(1:h,1);
L = l*opt.L;
opts.isreal = false;
max_eigval = real(eigs(@Hfun, h, 1, 'LM', opts));
eigens = real(eigs(@Hfun_shifted, h, max_number, 'LM', opts)) + max_eigval;
% figure; plot(sort(eigens), '*');

    function v = Hfun(f_hat)
        G = [0:h/2-1 -h/2:-1]'*n/h + k;
        f = ifft(f_hat);
        v = -0.5*(1i*2*pi/L*G).^2.*f_hat + fft(V.*f);
    end
    function v = Hfun_shifted(f_hat)
        G = [0:h/2-1 -h/2:-1]'*n/h + k;
        f = ifft(f_hat);
        v = (-0.5*(1i*2*pi/L*G).^2 - max_eigval).*f_hat + fft(V.*f);
    end
%     function v = afun(f_hat)
%         G = [0:h/2-1 -h/2:-1]' + k;
%         f = ifft(f_hat);
%         v = -0.5*(1i*2*pi/L*G).^2.*f_hat + fft(V.*f);
%     end
%     function v = bfun(f_hat)
%         G = [0:h/2-1 -h/2:-1]' + k;
%         f = ifft(f_hat);
%         v = (-0.5*(1i*2*pi/L*G).^2 - max_eigval).*f_hat + fft(V.*f);
%     end

end