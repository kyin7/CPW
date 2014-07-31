function f2 = spectral_diff(L, N, f, n)
% y = linspace(0,L,N+1); x = y(1:N);
k = [0:N/2-1 -N/2:-1]';
if mod(n,2)==1
    k(N/2+1) = 0;
end
f1 = fft(f);
f2 = real(ifft(((1i*2*pi/L*k).^n.*f1))); % f2 = d^n f
end