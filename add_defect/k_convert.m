function f2 = k_convert(f, N, k)
% f is an array of complex numbers representing Fourier coefficients of
% some real one variable function F at  wave numbers k + G, where G is the 
% 1-D reciprocal lattice N/h*(-h/2:h/2-1), and k can be 0, 1, ..., N/(2*h)
% when k = 0 or N/(2*h), k + G is symmetric aboud zero, so f shoud be skew-
% symmetric, i.e. f(omega) = conj(f(-omega)) for omega in k + G. In these
% two cases, we covert f by separating the real and imaginary part of it.
h = length(f);
if k==0
    f2 = [f(1);...
          sqrt(0.5)*(f(2:h/2) + f(h:-1:h/2+2));...
          f(h/2+1);...
          sqrt(0.5)*1i*(-f(2:h/2)+f(h:-1:h/2+2))];
elseif k==N/(2*h)
    f2 = sqrt(0.5)*[f(1:h/2) + f(h:-1:h/2+1); ...
                  1i*(-f(1:h/2)+f(h:-1:h/2+1))];
else
    f2 = f;
end
end