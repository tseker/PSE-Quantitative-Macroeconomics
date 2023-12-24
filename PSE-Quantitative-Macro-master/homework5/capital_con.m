function [c,ceq] = capital_con(x,gamma,theta, delta ,phi,beta, alpha, P, Z, kgrid,k_ind,z_ind)
    kprime = exp(Z(z_ind))*kgrid(k_ind)^alpha*x(2)^(1-alpha)+(1-delta)*kgrid(k_ind)-x(1);
    c(1) = kprime - max(kgrid);
    c(2) = min(kgrid)-kprime ;
    ceq =[];
end