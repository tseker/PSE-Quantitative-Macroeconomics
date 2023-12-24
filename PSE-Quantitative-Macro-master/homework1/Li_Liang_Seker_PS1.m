% ==========================
% Problem set solving the stochastic neoclassical growth model using a
% sequence method
% ==========================
% clear the workspace
clear
% close all figures
close all
cd 'D:\PSE\M2S1\quantative macro\homework\homework1'
addpath('D:\PSE\M2S1\quantative macro\homework\homework1')


%% pre setting
% ============
% parameters  - you may have to change this according to instructions
% ============
alpha=0.4; % capital share
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 1.0001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.025;

% ============
% options and convergence criteria
% ============

criter_V = 1e-7; % conv criterion for value function
T=150; % periods for transition

%mean of capital non-stochastic steady state
kbar = (alpha*beta/(1-beta*(1-delta)))^(1/(1-alpha));
% initial level of capital in the transition
k_0= 0.9*kbar; % you may have to change this from the problem set instructions

%steady-state of consumption
cbar = kbar^alpha - delta * kbar;



%% broyden
% ==============
% Solve deterministic sequence
% ==============

% as usual we need an initial guess for the sequences of k and c - here we
% just use steady state
x0=[kbar*ones(T,1)*0.9;cbar*ones(T,1)*0.9];
x0(1) = k_0;
x0(2*T) = cbar;


zpath=ones(T,1);
zpath(1)=1.01; % log AR(1) Process
for i=2:T-1
zpath(i)=exp(log(zpath(i-1)));
end

% now we need a function that returns errors of the equation system,
% constraining k(0)=k_0, and c(T)=cbar

% ==============
% a. Broyden's method
% ==============

%Initial guess for jacobian - use finite difference at steady state
%sequences of k and c

clear dx J
J = eye(2*T,2*T);
for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    J(:,i)=(rbc_obj_start_updated(x0+dx,alpha,beta,sigma,delta,k_0,cbar,ones(T,1))-...
            rbc_obj_start_updated(x0,alpha,beta,sigma,delta,k_0,cbar,ones(T,1)))./dx(i);
end

crit=1e-5;
x=x0;
f=rbc_obj_start_updated(x,alpha,beta,sigma,delta,k_0,cbar,zpath);

while max(abs(f))>crit

dx =  - pinv(J)*rbc_obj_start_updated(x,alpha,beta,sigma,delta,k_0,cbar,zpath);
x=x+dx; % x_n+1
f = rbc_obj_start_updated(x,alpha,beta,sigma,delta,k_0,cbar,ones(T,1));
J = J + (f * dx')/(dx'*dx) ;
end

k_trans_br=x(1:T,1);
c_trans_br=x(T+1:2*T,1);

%or we just let matlab solve it

[x,fsol]=fsolve(@(x)rbc_obj_start_updated(x,alpha,beta,sigma,delta,k_0,cbar,ones(T,1)),x0 );

k_trans=x(1:T,1);
c_trans=x(T+1:2*T,1);


%% multiple shooting
% ==============
% b. Multiple shooting using bisection
% ==============
c_0l=cbar*0.9*k_0/kbar;
c_0u=cbar*1.8*k_0/kbar;

% now you need to write a function that caculates cT using the dynamic
% equations and starting at guess c0, given k0
% 
% note that the equations may
% produce errors on the way if c or k go to 0
% you may have to redefine the errors to accomodate this
errorl=cT(c_0l,k_0,T,alpha,delta,sigma,beta)-cbar; 
erroru=cT(c_0u,k_0,T,alpha,delta,sigma,beta)-cbar;
while abs(c_0l-c_0u)>criter_V
    if errorl>0 || erroru<0
        disp('bisection not working')
        break
    end

    c_0n=(c_0l+c_0u)/2;
    errorn=cT(c_0n,k_0,T,alpha,delta,sigma,beta)-cbar;
    if errorn>0
        c_0u=c_0n;
    else
        c_0l=c_0n;
    end
end

c0 = c_0n;



% c0 =  cbar*1.8*k_0/kbar;
% 
% c0 = fzero(@(x)(cT(x,k_0,T,alpha,delta, sigma, beta)-cbar), c0 );



%% find the path

c_path_ms = zeros(T,1);
k_path_ms = zeros(T,1);
k_path_ms(1) = k_0;
c_path_ms(1)= c0;
    c = c0;
    k = k_0;
    
    for i = 1:T-1
       k_p =k;
       c_p = c;
       k = k_p^alpha+(1-delta) * k_p - c_p;
       if k <=0 
           k = k_p;
           break
       end
       c = ( 1/beta* c_p^(-sigma)/ (alpha * k^(alpha-1) +(1-delta)))^(-1/sigma);
       if ~isreal(c)
           c = c_p;
          break 
       end
       k_path_ms(i+1) = k;
       c_path_ms(i+1) = c;
    end






%% finitely living agents

% ==============
% 4. finetly living agents
% ==============


% k_{T} is zero, according to the transversality condition
% So for the Euler Equation at Perod T-1, we have 
% u'(k_{T-1}^alpha + (1-delta)*k_{T-1}-k_{T}) = beta* u'(k_T^alpha +(1-delta)*k_T) *(1-delta)
% We can either use the reverse shooting(guess and shoot from k_T to match k0)
% Or we can use the Broyden's method, changing a little bit of the equation
% Here, we choose the more convenient way, which is the Broyden's method.

% Use an adjusted rbc_obj function to accomodate finite living problem.

% ==============
% 4. finetly living agents
% ==============


% k_{T} is zero, according to the transversality condition
% So for the Euler Equation at Perod T-1, we have 
% u'(k_{T-1}^alpha * (1-delta)*k_T-k_{T+1}) = beta* u'(k_T^alpha + (1-delta)*k_T)
% We can either use the reverse shooting(guess and shoot from k_T to match k0)
% Or we can use the Broyden's method, changing a little bit of the equation
% Here, we choose the more convenient way, which is the Broyden's method.

% Use an adjusted rbc_obj function to accomodate finite living problem.

% ==============
% Using fsolve
% ==============

T = 10;
x0_T10 = [kbar*ones(T,1)*0.1;cbar*ones(T,1)*0.1];
x_finite_living_T10 = fsolve(@(x)(rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0)), x0_T10 );

k_finite_living_T10 = x_finite_living_T10(1:T);
c_finite_living_T10 = x_finite_living_T10(T+1:2*T);

T = 100;
x0_T100 = [kbar*ones(T,1)*0.1;cbar*ones(T,1)*0.1];
x_finite_living_T100 = fsolve(@(x)(rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0)), x0_T100 );

k_finite_living_T100 = x_finite_living_T100(1:T);
c_finite_living_T100 = x_finite_living_T100(T+1:2*T);

T = 200;
x0_T200 = [kbar*ones(T,1)*0.1;cbar*ones(T,1)*0.1];
x_finite_living_T200 = fsolve(@(x)(rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0)), x0_T200 );

k_finite_living_T200 = x_finite_living_T200(1:T);
c_finite_living_T200 = x_finite_living_T200(T+1:2*T);

%%
% ==============
% Using Broyden method
% ==============

% ========
% T = 10
% ========

clear dx J T x0 x

T = 10;
x0 = [kbar*ones(T,1)*0.9;cbar*ones(T,1)*0.9];
J = eye(2*T,2*T);

for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    J(:,i)=(rbc_obj_finite_living(x0+dx,alpha,beta,sigma,delta,k_0)-...
            rbc_obj_finite_living(x0,alpha,beta,sigma,delta,k_0))./dx(i);
end 


crit=1e-5;
x=x0;
f=rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);

while max(abs(f))>crit

dx =  - pinv(J)*rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);
x=x+dx; % x_n+1
f = rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);
J = J + (f * dx')/(dx'*dx) ;
end

k_trans_br_T10=x(1:T,1);
c_trans_br_T10=x(T+1:2*T,1);

% ========
% T = 100
% ========

clear dx J T x0 x

T = 100;
x0 = [kbar*ones(T,1)*0.9;cbar*ones(T,1)*0.9];
J = eye(2*T,2*T);

for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    J(:,i)=(rbc_obj_finite_living(x0+dx,alpha,beta,sigma,delta,k_0)-...
            rbc_obj_finite_living(x0,alpha,beta,sigma,delta,k_0))./dx(i);
end

crit=1e-5;
x=x0;
f=rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);

while max(abs(f))>crit

dx =  - pinv(J)*rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);
x=x+dx; % x_n+1
f = rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);
J = J + (f * dx')/(dx'*dx) ;
end

k_trans_br_T100=x(1:T,1);
c_trans_br_T100=x(T+1:2*T,1);

% ========
% T = 200
% ========

clear dx J T x0 x

T = 200;
x0 = [kbar*ones(T,1)*0.9;cbar*ones(T,1)*0.9];
J = eye(2*T,2*T);

for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    J(:,i)=(rbc_obj_finite_living(x0+dx,alpha,beta,sigma,delta,k_0)-...
            rbc_obj_finite_living(x0,alpha,beta,sigma,delta,k_0))./dx(i);
end

crit=1e-5;
x=x0;
f=rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);

while max(abs(f))>crit

dx =  - pinv(J)*rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);
x=x+dx; % x_n+1
f = rbc_obj_finite_living(x,alpha,beta,sigma,delta,k_0);
J = J + (f * dx')/(dx'*dx) ;
end

k_trans_br_T200=x(1:T,1);
c_trans_br_T200=x(T+1:2*T,1);


%% Reverse Multiple Shooting


%% multiple shooting, T = 10

kT_u=10;
kT_l=0;

% Reverse multiple Shooting! Match k0 with the initial value of the computed path
T = 10;
[k_path_u, c_path_u] = reverse_ms(kT_u, T,alpha,delta,sigma,beta);
[k_path_l, c_path_l] = reverse_ms(kT_l, T,alpha,delta,sigma,beta);
errorl = k_path_l(1)-k_0;
erroru = k_path_u(1)-k_0; 

while abs(kT_u-kT_l)>criter_V
    if errorl>0 || erroru<0
        disp('bisection not working')
        break
    end

    kT_n=(kT_l+kT_u)/2;
    [k_path_n, c_path_n] = reverse_ms(kT_n, T,alpha,delta,sigma,beta);
    errorn=k_path_n(1)-k_0;
    if errorn>0
        kT_u=kT_n;
    else
        kT_l=kT_n;
    end
end

k_path_ms_t10 = k_path_n;
c_path_ms_t10 = c_path_n;







%% multiple shooting T = 100

kT_u=7;
kT_l=0;

% Reverse multiple Shooting! Match k0 with the initial value of the computed path
T = 100;
[k_path_u, c_path_u] = reverse_ms(kT_u, T,alpha,delta,sigma,beta);
[k_path_l, c_path_l] = reverse_ms(kT_l, T,alpha,delta,sigma,beta);
errorl = k_path_l(1)-k_0;
erroru = k_path_u(1)-k_0; 

while abs(kT_u-kT_l)>criter_V
    if errorl>0 || erroru<0
        disp('bisection not working')
        break
    end

    kT_n=(kT_l+kT_u)/2;
    [k_path_n, c_path_n] = reverse_ms(kT_n, T,alpha,delta,sigma,beta);
    errorn=k_path_n(1)-k_0;
    if errorn>0
        kT_u=kT_n;
    else
        kT_l=kT_n;
    end
end

k_path_ms_t100 = k_path_n;
c_path_ms_t100 = c_path_n;


% c0 =  cbar*1.8*k_0/kbar;
% 
% c0 = fzero(@(x)(cT(x,k_0,T,alpha,delta, sigma, beta)-cbar), c0 );


%% multiple shooting, T = 200


kT_u=7;
kT_l=0;

% Reverse multiple Shooting! Match k0 with the initial value of the computed path
T = 200;
[k_path_u, c_path_u] = reverse_ms(kT_u, T,alpha,delta,sigma,beta);
[k_path_l, c_path_l] = reverse_ms(kT_l, T,alpha,delta,sigma,beta);
errorl = k_path_l(1)-k_0;
erroru = k_path_u(1)-k_0; 

while abs(kT_u-kT_l)>criter_V
    if errorl>0 || erroru<0
        disp('bisection not working')
        break
    end

    kT_n=(kT_l+kT_u)/2;
    [k_path_n, c_path_n] = reverse_ms(kT_n, T,alpha,delta,sigma,beta);
    errorn=k_path_n(1)-k_0;
    if errorn>0
        kT_u=kT_n;
    else
        kT_l=kT_n;
    end
end

k_path_ms_t200 = k_path_n;
c_path_ms_t200 = c_path_n;


%%
plotting
