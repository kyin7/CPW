function f2 = k_recover(f, N, k) 
% the inverse transform of k_convert
h = length(f);
if k==0
    f2 = [(f(1));...
         (f(2:h/2) + 1i*f(h/2+2:h))/sqrt(2);...
          (f(h/2+1));...
         (f(h/2:-1:2) - 1i*f(h:-1:h/2+2))/sqrt(2)];
elseif k==N/(2*h)
    f2 = [f(1:h/2) + 1i*f(h/2+1:h);...
          f(h/2:-1:1) - 1i*f(h:-1:h/2+1)]/sqrt(2);
else
    f2 = f;
end
end