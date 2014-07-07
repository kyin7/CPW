function Psi_hat_omega = Psi_hat_matrix(Psi_hat, omega, N, h, k)
p_shift = N/2 + 1;
if k==1
    Psi_hat_omega = zeros(h,1);
    return;
end
if omega==0 % size(Psi_hat_omega) = (floor(h/2)+1,1)
    p = N/h*(-floor(h/2-1/2) : -1);    
    if mod(h,2)==0
        Psi_hat_omega = [Psi_hat(-N/2+p_shift,1:k-1); sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            Psi_hat(p_shift,1:k-1); sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
    else
        Psi_hat_omega = [sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            Psi_hat(p_shift,1:k-1); sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
    end
elseif omega>=1 && omega<=floor(N/(2*h)-1/2) % size(Psi_hat_omega) = (floor(h/2),1)
    p = omega + (-floor(h/2):floor(h/2-1/2))*N/h;
    Psi_hat_omega = imag(Psi_hat(p + p_shift,1:k-1)*diag((1i).^(1:k-1)));   
elseif omega==N/(2*h)
    p = -N/h*floor(h/2) + N/(2*h) : N/h : -N/(2*h);
    if mod(h,2)==0
        Psi_hat_omega = [sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
    else
        Psi_hat_omega = [Psi_hat(-N/2+p_shift,1:k-1); sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
    end
end