function error=rbc_obj_start_updated(x,alpha,beta,sigma,delta,k0,cbar,z_trans)

%x(1:T) is capital
%x(T+1:2T) is consumption

T=length(x)/2;

error(1,1)=x(1,1)-k0;
error(2*T,1)=x(2*T,1)-cbar;

% calculate errors in budget constraint
for t=2:T
   % c_{t-1} - k_{t-1}^alpha -(1-delta)k_{t-1} + k_t
   error(t,1)=  x(T+t-1,1)-x(t-1,1)^alpha-(1-delta)*x(t-1,1)+x(t,1) ;
end

% calculate errors in EE
for t=1:T-1
      error(T+t,1)= ( x(T+t,1) )^(-sigma) -...
                    beta*( x(T+t+1,1)  )^(-sigma)*(alpha * x(t+1,1)^(alpha-1) + 1 - delta  );
end

end