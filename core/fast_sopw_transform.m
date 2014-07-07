function g_sopw = fast_sopw_transform(g_hat, N, h, k)
% N is even, h divides N
% size(g_sopw) = (N/h,1); % the coefficient vector for g under level-k SOPW
% g_hat = fftshift(fft(g));
omega_offset = N/2 + 1;
u_hat = zeros(N/h,1);
if floor(N/(2*h)) == N/(2*h)
    if k==1
        omega = -N/(2*h) : N/(2*h)-1;
        u_hat(1,1)= real(g_hat(omega(1)+omega_offset))/sqrt(h/2); 
        u_hat(2:N/h,1) = g_hat(omega(2:N/h)+omega_offset)/sqrt(h);
    elseif k==h
        omega = [-N/2 : -(h-1)*N/(2*h)-1 (h-1)*N/(2*h) : N/2-1];
        u_hat(1,1) = g_hat(omega(1)+omega_offset)/sqrt(h);
        u_hat(2:N/(2*h),1) = g_hat(omega(2:N/(2*h))+omega_offset)/((-1i)^(h-1)*sqrt(h));
        u_hat(N/(2*h)+1,1) = real(g_hat(omega(N/(2*h)+1)+omega_offset)/(1i^(h-1)*sqrt(h/2)));
        u_hat(N/(2*h)+2:N/h,1) = g_hat(omega(N/(2*h)+2:N/h)+omega_offset)/(1i^(h-1)*sqrt(h));
    else
        omega = [-k*N/(2*h) : -(k-1)*N/(2*h)-1 (k-1)*N/(2*h): k*N/(2*h)-1];
        u_hat(1,1) = real(g_hat(omega(1)+omega_offset)/((-1i)^(k-1)*sqrt(h/2)));
        u_hat(2:N/(2*h),1) = g_hat(omega(2:N/(2*h))+omega_offset)/((-1i)^(k-1)*sqrt(h));
        u_hat(N/(2*h)+1,1) = real(g_hat(omega(N/(2*h)+1)+omega_offset)/(1i^(k-1)*sqrt(h/2)));
        u_hat(N/(2*h)+2:N/h,1) = g_hat(omega(N/(2*h)+2:N/h)+omega_offset)/(1i^(k-1)*sqrt(h));
    end
    if mod(k,2)==0
        g_sopw = ifft(u_hat);
    else
        g_sopw = (ifft(fftshift(u_hat)));
    end
else % floor(N/(2*h)) < N/(2*h) 
    if k==1
        omega = [0 : N/(2*h)-1/2 -N/(2*h)+1/2 : -1];
        u_hat = g_hat(omega+omega_offset)/sqrt(h);
    elseif k==h
        omega = [-N/2: -((h-1)*N/(2*h)+1/2) (h-1)*N/(2*h)+1/2 : N/2-1];
        u_hat(1,1) = g_hat(omega(1)+omega_offset)/sqrt(h);
        u_hat(2:N/(2*h)+1/2,1) = g_hat(omega(2:N/(2*h)+1/2)+omega_offset)/((-1i)^(h-1)*sqrt(h));     
        u_hat(N/(2*h)+3/2:N/h) = g_hat(omega(N/(2*h)+3/2:N/h)+omega_offset)/(1i^(h-1)*sqrt(h));
    elseif mod(k,2)==1 % k is odd
        omega = [(k-1)*N/(2*h): k*N/(2*h)-1/2 -k*N/(2*h)+1/2 : -(k-1)*N/(2*h)-1];
        u_hat(1,1) = real(g_hat(omega(1)+omega_offset)/(1i^(k-1)*sqrt(h/2)));
        u_hat(2:N/(2*h)+1/2,1) = g_hat(omega(2:N/(2*h)+1/2)+omega_offset)/(1i^(k-1)*sqrt(h));
        u_hat(N/(2*h)+3/2:N/h,1) = g_hat(omega(N/(2*h)+3/2:N/h)+omega_offset)/((-1i)^(k-1)*sqrt(h));
    else
        omega = [-k*N/(2*h) : -(k-1)*N/(2*h)-1/2 (k-1)*N/(2*h)+1/2 : k*N/(2*h)-1];
        u_hat(1,1) = real(g_hat(omega(1)+omega_offset)/((-1i)^(k-1)*sqrt(h/2)));
        u_hat(2:N/(2*h)+1/2,1) = g_hat(omega(2:N/(2*h)+1/2)+omega_offset)/((-1i)^(k-1)*sqrt(h));
        u_hat(N/(2*h)+3/2:N/h,1) = g_hat(omega(N/(2*h)+3/2:N/h)+omega_offset)/(1i^(k-1)*sqrt(h));
    end 
    if mod(k,2)==1
        g_sopw = (ifft(u_hat));
    else
        g_sopw = ifft(u_hat);
    end
end
% g_hat_omega = g_hat(omega + omega_shift)';
% if k==1
%     for a=0:N/h-1
%         g_sopw(a+1,1) = dot(g_hat_omega(end), u_hat(end,1)) +...
%             2*real(dot(g_hat_omega(1:end-1), u_hat(1:end-1,1).*exp(1i*2*pi*omega(1:end-1)'*a*h/N)));
%     end
% elseif k==h
%     for a=0:N/h-1
%         g_sopw(a+1,1) = dot(g_hat_omega(1), u_hat(1,1)) +...
%             2*real(dot(g_hat_omega(2:end), u_hat(2:end,1).*exp(1i*2*pi*omega(2:end)'*a*h/N)));
%     end
% else
%     for a=0:N/h-1
%         g_sopw(a+1,1) = 2*real(dot(g_hat_omega, u_hat.*exp(1i*2*pi*omega'*a*h/N)));
%     end
% end
% g_sopw = g_sopw*h/N;