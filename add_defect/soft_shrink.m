function u = soft_shrink(v, sigma)
u = max((abs(v) - sigma)./abs(v), 0) .* v;
