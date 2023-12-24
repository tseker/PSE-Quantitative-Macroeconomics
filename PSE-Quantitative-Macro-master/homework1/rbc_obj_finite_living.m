function error=rbc_obj_finite_living(x,alpha,beta,sigma,delta,k0)

%x(1:T) is capital, k0, k1, ..., kT-1, kT = 0, because of transversaltiy condition

T=length(x)/2;

% here, we don't have cbar or kbar in the steady state. So only need an initial capital
error(2*T,1)=x(1,1)-k0;


% calculate errors in budget constraint
for t=1:T
   % c_{t-1} - k_{t-1}^alpha -(1-delta)k_{t-1} + k_t
   if t == T
        error(t,1) = x(T+t,1)-x(t,1)^alpha-(1-delta)*x(t,1); % x(T+1), which is kT, equals 0
   else
        error(t,1)=  x(T+t,1)-x(t,1)^alpha-(1-delta)*x(t,1)+x(t+1,1) ;
   end
end

% calculate errors in EE
for t=1:T-1
   error(T+t,1)= ( x(T+t,1) )^(-sigma) -...
                   beta*( x(T+t+1,1)  )^(-sigma)*(alpha * x(t+1,1)^(alpha-1) + 1 - delta  );
end

end