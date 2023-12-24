function [terminal_val_c, terminal_val_k] = cT(c0, k0,T,alpha,delta,sigma,beta)
    c = c0;
    k = k0;
    for i = 1:T
       k_p =k;
       c_p = c;
       k = k_p^alpha+(1-delta) * k_p - c_p;
       if k <=0 
           break
       end
       c = ( 1/beta* c_p^(-sigma)/ (alpha * k^(alpha-1) +(1-delta)))^(-1/sigma);
       if ~isreal(c)
          break 
       end
    end
    terminal_val_c = c;
    terminal_val_k = k;
end


asdfasdfasdf


asdfasdfasdfasdfasdf