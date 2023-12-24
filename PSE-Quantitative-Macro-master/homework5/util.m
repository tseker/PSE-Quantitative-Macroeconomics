function newV = util(x,V ,gamma,theta, delta ,phi,beta, alpha, P, Z, kgrid ,k_ind, z_ind)

    kprime = exp(Z(z_ind))*kgrid(k_ind)^alpha*x(2)^(1-alpha)+(1-delta)*kgrid(k_ind)-x(1);
    [xx,yy] = ndgrid(real(kprime), exp(Z));

    newV =  (x(1)^(1-gamma)-1)/(1-gamma) - theta/(1+phi)*x(2)^(1+phi)...
               + beta* (P(z_ind,:)*V(xx,yy)');
 end