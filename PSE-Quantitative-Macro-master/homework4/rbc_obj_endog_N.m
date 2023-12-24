function   error =  rbc_obj_endog_N(x,alpha,beta,gamma,delta,theta,psi,k_0,cbar,kbar,Lbar,z_trans)
    % x: inputs , capital,  consumption path, labor path, T* 3
    % outputs: a 3T times 1 vector with errors of the Euler equation for
    % capital investment, intratemporal condition for labor supply
    % and the feasibility constraint
    
    %x(1:T) is capital, k0, k1, ..., kT-1, put kT = kbar, because of transversaltiy condition
    
    s = size(x);
    
    T=s(1);
    
    error = zeros(3*T,1); % 1 to T, budget constraint; T to 2T,
    error(T)=x(T,1)-kbar;
    error(2*T)=x(1,1)-k_0;
    
% calculate errors in budget constraint
for t=1:T-1
    % ct + Kt+1 - A_t Kt^alpha N_t^{1-alpha} - (1 - delta)*Kt = 0;
   error(t)=  x(t,2)+ x(t+1,1) - z_trans(t)*x(t,1)^(alpha)* x(t,3)^(1-alpha) - (1-delta)* x(t,1) ;
end

% calculate errors in EE
for t=1:T-1
    % Ct^(-gamma) - beta * Ct+1^(-gamma) * ( alpha * At+1 * Kt+1^(alpha-1)* Nt+1^(1-alpha)+ 1 - delta )
   error(T+t)= x(t,2)^(-gamma) - beta * x(t+1,2)^(-gamma) * ( alpha * z_trans(t+1) * x(t+1,1)^(alpha-1)* x(t+1,3)^(1-alpha)+ 1 - delta )  ;
end
% errors in FOC for labor supply
for t=1:T
    %theta Nt^psi - Ct^{-gamma} (1-alpha) At Kt^alpha Nt^{-alpha}
   error(2*T+t)= theta * x(t,3)^(psi) - x(t,2)^(-gamma) * (1-alpha) * z_trans(t)* x(t,1)^(alpha)*x(t,3)^(-alpha) ;
end

end