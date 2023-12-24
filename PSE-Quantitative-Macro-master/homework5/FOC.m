function error=FOC(x,K,Kprime,z,gamma,theta,phi,alpha,beta,delta)
% This is the error of budget constraint and the error of the
% intra-temporal labor-consumption substitution
% x is the consumption and labor

c = x(1);
L = x(2);
error = zeros(2,1);
error(1) = c + Kprime - z*K^alpha*L^(1-alpha) - (1-delta)* K;
error(2) = theta*L^(phi+alpha)  - c^(-gamma)*(1-alpha)*z*K^alpha ;


end