function test_conv(opt)
% n = 16;
% u = randn(n,1);
% v = randn(n,1);
% z = convolve(u,v);
% u_hat = fft(u);
% v_hat = fft(v);
% w_hat = convolve(u_hat,v_hat);
% w = w_hat/n - fft(u.*v);

V = opt.V;
V_hat = fft(V);
N = opt.N;
h = opt.h;
G = (0:h-1)*N/h;
V_bar = V_hat(G+1);
f = randn(N,1);
u = V.*f;
v_hat = zeros(N,1);
f_hat = fft(f);
for k=0 : N/(2*h)
    f_k = f_hat(k + G + 1);
    v_k = convolve(V_bar, f_k)/N;
    v_hat(k + G + 1) = v_k;
end
for k= N/(2*h)+1 : N/h-1
    v_hat(k + G + 1) = conj(v_hat(-k - G + N + 1));
end
v = ifft((v_hat));
w = u - v;
w_hat = fft(u) - v_hat;
end