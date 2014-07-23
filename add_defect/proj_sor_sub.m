function x_p = proj_sor_sub(v, k, opt, omega, P_omega)
% x_p = argmin (x-v).^2 s.t. x is in Range(P_omega) and
% ||x||^2=h
% n = size(P_omega,1);
N = opt.N;
h = opt.h;
p_shift = N/2+1;
if omega==0
    p = N/h*(-floor(h/2-1/2) : -1);
    if mod(h,2)==0
        vp = [v(-N/2+p_shift,1); sqrt(2)*real(v(p + p_shift,1));...
            v(p_shift,1); sqrt(2)*imag(v(p + p_shift,1))];
        Dp = [1; sqrt(2)*ones(floor(h/2-1/2),1); 1; sqrt(2)*ones(floor(h/2-1/2),1)];
    else
        vp = [sqrt(2)*real(v(p + p_shift,1));...
            v(p_shift,1); sqrt(2)*imag(v(p + p_shift,1))];
        Dp = [sqrt(2)*ones(floor(h/2-1/2),1); 1; sqrt(2)*ones(floor(h/2-1/2),1)];
    end
    
elseif omega>=1 && omega<=floor(N/(2*h)-1/2)
    p = (omega - floor(h/2)*N/h) : N/h : (omega + floor(h/2-1/2)*N/h);
    if mod(k,2) == 0
        vp = v(p + p_shift,1);
%         vp = imag(v(p + p_shift,1));
    else
        vp = v(p + p_shift,1);
%         vp = real(v(p + p_shift,1));
    end
    Dp = ones(h,1);
elseif omega==N/(2*h)
    p = -N/h*floor(h/2) + N/(2*h) : N/h : -N/(2*h);    
    if mod(h,2)==1
        vp = [v(-N/2+p_shift,1); sqrt(2)*real(v(p + p_shift,1));...
            sqrt(2)*imag(v(p + p_shift,1))];
        Dp = [1; sqrt(2)*ones(floor(h/2),1); sqrt(2)*ones(floor(h/2),1)];
    else
        vp = [sqrt(2)*real(v(p + p_shift,1));...
            sqrt(2)*imag(v(p + p_shift,1))];
        Dp = [sqrt(2)*ones(floor(h/2),1); sqrt(2)*ones(floor(h/2),1)];
    end
end
% yp = P_omega*vp;
yp = vp - P_omega * (P_omega'*vp)/h; % G-S orthogonalization to previous level modes
if k<0*opt.max_level
if norm(yp)>1e-3*sqrt(h)
    yp = yp/norm(yp)*sqrt(h);
else
    yp = 0*yp;
end
end
yp = yp./Dp;
if omega==0
    if mod(h,2)==0
        x_p = [yp(1,1); yp(2:floor(h/2-1/2)+1,1) + 1i*yp(floor(h/2-1/2)+3:h,1);...
            yp(floor(h/2-1/2)+2,1)];
%         if mod(k,2)==1
%             x_p = [yp(1,1); yp(2:floor(h/2-1/2)+1,1); yp(floor(h/2-1/2)+2,1)];
%         else
%             x_p = [yp(1,1); 1i*yp(floor(h/2-1/2)+3:h,1); yp(floor(h/2-1/2)+2,1)];
%         end
    else
        x_p = [yp(1:floor(h/2-1/2),1) + 1i*yp(floor(h/2-1/2)+2:h,1); yp(floor(h/2-1/2)+1,1)];
%         if mod(k,2)==1
%             x_p = [yp(1:floor(h/2-1/2),1); yp(floor(h/2-1/2)+1,1)];
%         else
%             x_p = [1i*yp(floor(h/2-1/2)+2:h,1); yp(floor(h/2-1/2)+1,1)];
%         end
    end
elseif omega>=1 && omega<=floor(N/(2*h)-1/2)
    if mod(k,2)==0
        x_p = yp;
%         x_p = 1i*yp;
    else
        x_p = yp;
    end
elseif omega==N/(2*h)
    if mod(h,2)==1
        x_p = [yp(1,1); yp(2:floor(h/2)+1,1) + 1i*yp(floor(h/2)+2:h,1)];
    else
        x_p = yp(1:floor(h/2),1) + 1i*yp(floor(h/2)+1:h,1);
    end
end
        
    
    
