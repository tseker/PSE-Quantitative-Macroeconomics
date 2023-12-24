function [k_path, c_path] = reverse_ms(kT, T,alpha,delta,sigma,beta)
    % kT+1 = 0
    k_path = zeros(T,1);
    c_path = zeros(T,1);
    c_path(T) = kT^alpha+(1-delta)*kT;
    k_path(T) = kT;
    
    for t = flip(1:T-1)
        c_path(t) = (beta* c_path(t+1)^(-sigma) * (1-delta + alpha * k_path(t+1)^(alpha-1)))^(-1/sigma);
        k_path(t) = fzero(@(x)( x^alpha + (1-delta)*x - k_path(t+1)- c_path(t)),[0,2*k_path(t+1)]  );
        if k_path(t) <1e-3 
            break
        end
    end
end

