function v = L1_step(w, opt)
if opt.mu == 0
    v = w;
else
    mu = opt.mu*sqrt(opt.N/opt.L);
    v = soft_shrink(w, 1/(mu*opt.gamma));
end