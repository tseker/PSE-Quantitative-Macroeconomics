function [policy_fun_c, policy_fun_kprime, value_fun] = ValueFunIteration_Discrete( N, theta, delta, beta, sigma, criter_V, Howard )
    
   kbar =((1/beta-1+delta)/(theta))^(1/(theta-1));
    
    linear = 1;
    if delta==1
        if linear==1
            kgrid=linspace(kbar/2,2*kbar,N); %kgrid vector is created using linspace to generate N equally spaced values between kbar/2 and 2*kbar
        else
            temp=linspace(0,0.5,N).^5/0.5^5*(2*kbar-kbar/2); %This operation results in a non-linear grid where values are closer together near the lower end of the interval (0 to 0.5) and spaced further apart as you approach 0.5.
            kgrid=kbar/2+temp; %non-linear grid of values that starts at kbar/2 and becomes increasingly spaced as it moves toward 2*kbar
        end
        else
        kgrid=linspace(3*kbar/4 ,5*kbar/4,N); %linear grid of values starting at kbar/4 and ending at 2*kbar
    end
    
    
    
    
    
    
    % Value function initialization
    
    c = zeros(N,N);
    u = zeros(N,N);
    
    if sigma ~=1
        for i=1:N
            for j=1:N
                c(i,j)= kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j); %Feasibility constraint
                if c(i,j)>0 && c(i,j) <= kgrid(i)^theta
                    u(i,j)=(c(i,j)^(1-sigma)-1)/(1-sigma); %Utility function
                else
                    u(i,j)=-1e50*((kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j))<=0);
                end
            end
        end
    else
        for i=1:N
            for j=1:N
                c(i,j)= kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j); %Feasibility constraint
                if c(i,j)>0 && c(i,j) <= kgrid(i)^theta
                    u(i,j)=log(c(i,j)); %Utility function
                else
                    u(i,j)=-1e50*((kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j))<=0);
                end
            end
        end
    end
    
    
    dV=1; % control the convergence criterion in the value function iteration (?)
    V=zeros(N,1);
    iter=0;
    tic

    % v. main loop
    VV = zeros(N,N);
    Vnew = zeros(N,1);
    ind = zeros(N,1);
    kprime_vf = zeros(N,1);
    c_vf = zeros(N,1);
    
    while dV>criter_V
        iter=iter+1;
        for i=1:N % loop over capital today
            for j=1:N %... and over capital tomorrow
                if kgrid(j) >= (1-delta) * kgrid(i)  && kgrid(j) <= (1-delta)*kgrid(i)+kgrid(i)^theta
                    VV(i,j)=u(i,j)+beta*V(j); % calculate value for each i,j
                else
                    VV(i,j) = -1e30;
                end
            end
            % take maximum over capital tomorrow
            [Vnew(i),ind(i)] = max(VV(i,:));
            % record policy function
            kprime_vf(i) = kgrid(ind(i));
            c_vf(i) = c(i,ind(i));
        end
        % Howard - doesn't help much here
        % In standard policy iteration, the policy is updated for all states in every iteration. 
        % In Howard's method, only states where the policy can be improved are updated. 
        % This selective update can significantly reduce the number of iterations required to find the optimal policy.
        if Howard==1 && iter>3
            dVV=1;
            while dVV>criter_V
                for i=1:N
                    % v(k) = u(k)+beta*v(k')
                    Vnewhow = (c_vf.^(1-sigma)-1)/(1-sigma) + beta*Vnew(ind);
                    clear temp
                end
                dVV=max(max(abs(Vnewhow-Vnew)));
                %disp(dVV)
                Vnew=Vnewhow;
            end
        end

        % calculate convergence criterion
        dV=max(abs(Vnew-V));

        % updated value function and one time return fuction u
        V=Vnew;
%         disp('dV')
%         disp(dV)
    end
    
    policy_fun_c = c_vf;
    policy_fun_kprime = kprime_vf;
    value_fun=V;
    
end