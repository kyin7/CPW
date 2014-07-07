function x_p = solve_x(v, lambda, r, N, h, k, omega, P_omega)
% x_p = argmin a.*x.^2 + lambda*(x-v).^2 s.t. x is in Range(P_omega) and
% ||x||^2=h
% n = size(P_omega,1);
p_shift = N/2+1;
if omega==0
    p = N/h*(-floor(h/2-1/2) : -1);
    if mod(h,2)==0
        vp = [v(-N/2+p_shift,1); sqrt(2)*real(v(p + p_shift,1));...
            v(p_shift,1); sqrt(2)*imag(v(p + p_shift,1))];
        ap = [-N/2 p 0 p].^2;
        Dp = [1; sqrt(2)*ones(floor(h/2-1/2),1); 1; sqrt(2)*ones(floor(h/2-1/2),1)];
    else
        vp = [sqrt(2)*real(v(p + p_shift,1));...
            v(p_shift,1); sqrt(2)*imag(v(p + p_shift,1))];
        ap = [p 0 p].^2;
        Dp = [sqrt(2)*ones(floor(h/2-1/2),1); 1; sqrt(2)*ones(floor(h/2-1/2),1)];
    end
    
elseif omega>=1 && omega<=floor(N/(2*h)-1/2)
    p = (omega - floor(h/2)*N/h) : N/h : (omega + floor(h/2-1/2)*N/h);
    if mod(k,2) == 0
        vp = imag(v(p + p_shift,1));
    else
        vp = real(v(p + p_shift,1));
    end
    ap = p.^2;
    Dp = ones(h,1);
elseif omega==N/(2*h)
    p = -N/h*floor(h/2) + N/(2*h) : N/h : -N/(2*h);    
    if mod(h,2)==1
        vp = [v(-N/2+p_shift,1); sqrt(2)*real(v(p + p_shift,1));...
            sqrt(2)*imag(v(p + p_shift,1))];
        ap = [-N/2 p p].^2;
        Dp = [1; sqrt(2)*ones(floor(h/2),1); sqrt(2)*ones(floor(h/2),1)];
    else
        vp = [sqrt(2)*real(v(p + p_shift,1));...
            sqrt(2)*imag(v(p + p_shift,1))];
        ap = [p p].^2;
        Dp = [sqrt(2)*ones(floor(h/2),1); sqrt(2)*ones(floor(h/2),1)];
    end
end
% xp = vp;
yp = zeros(h,1);
z = zeros(h,1);
jj_max = 10;
% epsilon0 = 1e-6;
% r = 1e5;
numer = lambda*vp;
denom = ap'+lambda+r;
for jj=1:jj_max 
    xp = (numer + r*(yp-z))./denom;
    yp = P_omega*(xp+z);
    yp = yp/norm(yp)*sqrt(h);
    z = z + xp - yp;
%     if norm(xp - yp)/norm(xp)<epsilon0 
%         break;
%     end
end
yp = yp./Dp;
if omega==0
    if mod(h,2)==0
%         x_p = [yp(1,1); yp(2:floor(h/2-1/2)+1,1) + 1i*yp(floor(h/2-1/2)+3:h,1);...
%             yp(floor(h/2-1/2)+2,1)];
        if mod(k,2)==1
            x_p = [yp(1,1); yp(2:floor(h/2-1/2)+1,1); yp(floor(h/2-1/2)+2,1)];
        else
            x_p = [yp(1,1); 1i*yp(floor(h/2-1/2)+3:h,1); yp(floor(h/2-1/2)+2,1)];
        end
    else
%         x_p = [yp(1:floor(h/2-1/2),1) + 1i*yp(floor(h/2-1/2)+2:h,1);...
%             yp(floor(h/2-1/2)+1,1)];
        if mod(k,2)==1
            x_p = [yp(1:floor(h/2-1/2),1); yp(floor(h/2-1/2)+1,1)];
        else
            x_p = [1i*yp(floor(h/2-1/2)+2:h,1); yp(floor(h/2-1/2)+1,1)];
        end
    end
elseif omega>=1 && omega<=floor(N/(2*h)-1/2)
    if mod(k,2)==0
        x_p = 1i*yp;
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
        
    
    
% if omega==0
%     x_p = real(xp);
% end
% vp = zeros(h,1);
% ap = zeros(h,1);
% if omega==0
%     real_index = 1:floor(h/2)+1;
%     imag_index = floor(h/2)+2:h;
%     cplx_index = 1+floor(h/2)-floor(h/2-1/2):floor(h/2);
%     vp(real_index,1) = real(v_p);
%     vp(imag_index,1) = imag(v_p(cplx_index));
%     ap(real_index,1) = a_p;
%     ap(imag_index,1) = a_p(cplx_index);
% elseif omega==N/(2*h)
%     real_index = 1:floor(h/2+1/2);
%     imag_index = floor(h/2+1/2)+1:h;
%     cplx_index = floor(h/2+1/2)-floor(h/2)+1:floor(h/2+1/2);
%     vp(real_index,1) = real(v_p);
%     vp(imag_index,1) = imag(v_p(cplx_index));
%     ap(real_index,1) = a_p;
%     ap(imag_index,1) = a_p(cplx_index);
% else
%     vp = v_p;
%     ap = a_p;
% end
% xp = vp;
% x_p = zeros(h,1);
% z = zeros(h,1);
% jj = 0; 
% jj_max = 10;
% epsilon0 = 1e-6;
% % r = 1e5;
% while norm(xp - x_p)/norm(xp)>epsilon0 && jj<jj_max
%     xp = (lambda*vp + r*(x_p-z))./(ap+lambda+r);
%     x_p = P_omega*(xp+z);
%     x_p = x_p/norm(x_p)*sqrt(h);
%     z = z + xp - x_p;
%     jj = jj + 1;
% end
% xp = x_p;
% if omega == 0 || omega == N/(2*h)
%     x_p = xp(real_index,1);
%     x_p(cplx_index,1) = x_p(cplx_index,1) + 1i*xp(imag_index,1);
% else
%     x_p = xp;
% end
% end

