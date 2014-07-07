function [u, u_hat] = fourier_orthogonal(N, h, Psi_hat, k)
u_hat = zeros(N,1);
p_shift = N/2 + 1;
if k==1
    [u, u_hat] = sopwbasis(N, h, 1);
    u = ifftshift(u);
else
    for omega=0:floor(N/(2*h))
        if omega==0
            p = (-floor(k/2):-1)*N/h; 
            if k==h && mod(h,2)==0
                Psi_hat_omega = [Psi_hat(-N/2 + p_shift,1:k-1); sqrt(2)*real(Psi_hat(p(2:end) + p_shift,1:k-1));...
            Psi_hat(p_shift,1:k-1); sqrt(2)*imag(Psi_hat(p(2:end) + p_shift,1:k-1))];
                Q_hat_omega = null(Psi_hat_omega');
                u_hat(-N/2 + p_shift,1) = sqrt(h)*Q_hat_omega(1,1);
                u_hat(p(2:end) + p_shift,1) = sqrt(h/2)*(Q_hat_omega(2:h/2,1) +...
                    1i*Q_hat_omega(h/2+2:h,1));
                u_hat(p_shift,1) = sqrt(h)*Q_hat_omega(h/2+1,1);
                p_prime = (1:floor(h/2-1/2))*N/h;
                u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift,1));
            else
                Psi_hat_omega = [sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            Psi_hat(p_shift,1:k-1); sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
                Q_hat_omega = null(Psi_hat_omega');
                u_hat(p + p_shift,1) = sqrt(h/2)*(Q_hat_omega(1:floor(k/2),1) +...
                    1i*Q_hat_omega(floor(k/2)+2:2*floor(k/2)+1,1));
                u_hat(p_shift,1) = sqrt(h)*Q_hat_omega(floor(k/2)+1,1);
                p_prime = (1:floor(k/2))*N/h;
                u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift,1));
            end
            
        elseif omega>0 && omega<N/(2*h)
            p = omega + (-floor(h/2):floor(h/2-1/2))*N/h;
            [~, IX] = sort(p.^2);
            p_prime = p(IX(1:k));
            Q_hat_omega = null(Psi_hat(p_prime + p_shift,1:k-1)');
            u_hat(p_prime + p_shift,1) = sqrt(h)*Q_hat_omega;
            u_hat(-p_prime + p_shift, 1) = conj(u_hat(p_prime + p_shift,1));

        elseif omega==N/(2*h)
            p = (-floor(k/2-1/2):0)*N/h - N/(2*h);
            if k==h && mod(h,2)==1
                Psi_hat_omega = [Psi_hat(-N/2 + p_shift,1:k-1); sqrt(2)*real(Psi_hat(p(2:end) + p_shift,1:k-1));...
            sqrt(2)*imag(Psi_hat(p(2:end) + p_shift,1:k-1))];
                Q_hat_omega = null(Psi_hat_omega');
                u_hat(-N/2 + p_shift,1) = sqrt(h)*Q_hat_omega(1,1);
                u_hat(p(2:end) + p_shift,1) = sqrt(h/2)*(Q_hat_omega(2:h/2,1) +...
                    1i*Q_hat_omega(h/2+2:h,1));
                p_prime = N/(2*h) + N/h*(1:floor(h/2)-1);
                u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift,1));
            elseif k==h && mod(h,2)==0
                Psi_hat_omega = [sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
                Q_hat_omega = null(Psi_hat_omega');
                u_hat(p + p_shift,1) = sqrt(h/2)*(Q_hat_omega(1:floor(k/2+1/2),1) +...
                    1i*Q_hat_omega(floor(k/2+1/2)+1:2*floor(k/2+1/2),1));
                p_prime = -p;
                u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift,1));
            else
                Psi_hat_omega = [sqrt(2)*real(Psi_hat(p + p_shift,1:k-1));...
            sqrt(2)*imag(Psi_hat(p + p_shift,1:k-1))];
                Q_hat_omega = null(Psi_hat_omega');
                u_hat(p + p_shift,1) = sqrt(h/2)*(Q_hat_omega(1:floor(k/2+1/2),1) +...
                    1i*Q_hat_omega(floor(k/2+1/2)+1:2*floor(k/2+1/2),1));
                p_prime = (0:floor(k/2))*N/h + N/(2*h);
                u_hat(p_prime + p_shift,1) = conj(u_hat(-p_prime + p_shift,1));
            end

        end
    end
    u = ifft(fftshift(u_hat));
    u = ifftshift(u);
end
