function test_poisson
N = 100;
L = 4;
lambda = 1e2;
gamma = 1e2;
g = randn(N,1); g = g/norm(g);
u1 = Mfun1(g);
u2 = Mfun2(g);
err = norm(u1-u2)/norm(u1);

    function u = Mfun1(f) % u = (-1/2*Laplace + lambda + gamma)^(-1) f
        f1 = fft(f);
        k2 = [0:N/2-1 -N/2:-1]';
        D = (-0.5*(1i*2*pi/L*k2).^2+lambda+gamma);
        f2 = f1./D;
        u = real(ifft(f2));
    end
    function u = Mfun2(f) % u = (-1/2*Laplace + lambda + gamma)^(-1) f
        f1 = fft(f);
        dx = L/N;
        C = zeros(N,1);
        C(1,1) = 1/dx^2 + lambda + gamma;
        C(2,1) = -0.5/dx^2;
        C(N,1) = -0.5/dx^2;
        D = fft(C);
        f2 = f1./D;
        u = real(ifft(f2));
    end
end