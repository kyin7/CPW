function [v] = shift(u, j, h)
% v(x) = u(x+j*h)
N = length(u);
v = zeros(N,1);
jh = mod((1:N)+j*h-1,N);
v(1:N) = u(jh+1);
