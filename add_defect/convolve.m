function w = convolve(u,v)
% w is the convolution of u and v, where u, v are of the same length
% n = length(u);
u_hat = fft(u);
v_hat = fft(v);
w = ifft(u_hat.*v_hat);
end