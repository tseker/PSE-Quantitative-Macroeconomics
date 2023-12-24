
% ==========================
% Problem set solving the RBC model different
% ways
% ==========================
clear
% close all figures
close all
% change location to the folder where this m file is saved
mfile_name          = mfilename('fullpath');  % This line of code retrieves the full path of the current MATLAB file and assigns it to the variable mfile_name. The mfilename function with the argument 'fullpath' returns the full path of the currently executing file, including the file name and extension. The full path is then assigned to the variable mfile_name.
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

%% 
% ==============
% Calibration
% ==============

% ============
% Calibration targets
% ============

Kshare=1/3;
Rbar=1.01; % This is the interest rate minus depreciation
Lbar=1/3;
InvKrat=0.05; % depreciation, also the investment/capital ratio in the steady state

% ============
% Exogenous parameters
% ============
gamma = 2 ; % CRRA
rho = 0.95;   % persistence of TFP shock
sigmaepsilon = 0.007; % volatility of TFP shock
phi=1; %inverse Frisch elasticity


% Calibration
beta=1/Rbar;
delta=InvKrat;
Rbar = Rbar+delta; % Recover the true interest rate
alpha=Kshare;
margprod=1/beta-1+delta;
KNrat=(margprod/alpha)^(1/(alpha-1));
wbar=(1-alpha)*(margprod/alpha)^(alpha/(alpha-1));
cKrat = margprod/alpha-delta;
theta = wbar/(   Lbar^(gamma+phi)* ( (margprod/alpha) -delta)^gamma * (margprod/alpha)^(gamma/(alpha-1)) );
kbar = KNrat*Lbar;

% Steady state
cbar = kbar * cKrat;
ybar = kbar^alpha*Lbar^(1-alpha);


% ============
% options and convergence criteria
% ============
criter_V = 1e-3; % conv criterion for value function
N=50; % number of grid points
linear=1; % grid linear or not
M=5; % number of support points for the shock
T=1000; % periods for simulation
N_sim=200; % number of simulations



%% Discrete Markov Chain Simulation
% ==============
% Grids, transition probabilities, etc
% ==============
kmax=1.1*kbar;
kmin=kbar*0.9;
if delta==1
    if linear==1
        kgrid=linspace(kmin,kmax,N);
    else
        temp=linspace(0,0.5,N).^5/0.5^5*(kmax-kmin);
        kgrid=kmin+temp;
    end
else
    kgrid=linspace(kmin ,kmax,N);
end

% ============
% Markov chain
% ============
%[Z, P] = tauchen(M,0,rho,sigmaepsilon,1)
[Z,P] = rouwenhorst(M,0,rho,sigmaepsilon);

% simulate discrete Markov Process
rng(2);

z = zeros(N_sim, T);
for j=1:N_sim
    z(j,1)=5;
    for t=2:T
        z(j,t)=sum( cumsum(P(z(j,t-1),:))<rand(1)*ones(1,M)  )+1;
    end
end

% simulate continuous Markov Process

z_cont = zeros(N_sim, T);
for j=1:N_sim
    z_cont(j,1)=0;
    for t=2:T
        z_cont(j,t)=rho*z_cont(j,t-1)+normrnd(0,sigmaepsilon);
    end
end

Z_lev=exp(Z);
Z_sim=Z_lev(z); % Z_sim is the simulation for stochastic A (log(A)=z)


%% Value Function Iteration: Initialization
% ==============
% Option I:  Fully discretized value function iteration
% ==============

% In VFI, we need to use non-linear solver to solve for the choice
% variable, which is not very desirable.


options=optimset('MaxIter',1000,'MaxFunEval',1000);

% one period return
% c(k,z,k'), u(k,z,k'), L(k,z,k')
c = zeros(N,M,N);
L = zeros(N,M,N);
u = zeros(N,M,N);



tic
for k_ind = 1:N
    for z_ind = 1:M
        for kp_ind = 1:N
            if kgrid(kp_ind)<= exp(Z(z_ind))*kgrid(k_ind)^alpha * 1 + (1-delta)* kgrid(k_ind)
                % to ensure positive consumption when capital in the next period is too high
                temp= fsolve(@(x) (FOC(x,kgrid(k_ind),kgrid(kp_ind),exp(Z(z_ind)),gamma,theta,phi,alpha,beta,delta)  ),...
                               [cbar,Lbar], optimset('Display','off') ) ;
                c(k_ind,z_ind,kp_ind)= temp(1);
                L(k_ind,z_ind,kp_ind)= temp(2);
                u(k_ind,z_ind,kp_ind)= (c(k_ind,z_ind,kp_ind)^(1-gamma)-1)/(1-gamma) -...
                                        theta /(1+phi) * L(k_ind,z_ind,kp_ind)^(1+phi);
            else
                c(k_ind,z_ind,kp_ind) = -1; % force consumption to be negative when kprime is not feasible
                L(k_ind,z_ind,kp_ind) = 1;
                u = -Inf;
            end
        end
    end
end
toc


%% Value Function Iteration: discrete VFI
% iterate on value function
% initial guess
    uu = zeros(N,M);
    uu = uu-Inf;
    [x,y]=ndgrid(kgrid,exp(Z));
    cc= x.^alpha.*Lbar.^(1-alpha) - delta.*x;
    ll= (1/theta .* cc.^(-gamma) .* (1-alpha) .* x.^(alpha)).^(1/(phi+alpha));
    uu= (cc.^(1-gamma)-1)./(1-gamma) - theta./(1+phi) .* ll.^(1+phi);
    V = uu./(1-beta);

    

% iterate on Value Function
% V(k,z)
%iteration setting 
    iter=0;
    dV = 1;


tic
Vtemp = zeros(N,M,N);
Vnew = V;
c_pol_disc = zeros(N,M);
L_pol_disc = zeros(N,M);
kprime_pol_disc = zeros(N,M);
kprime_ind = zeros(N,M);

while dV>criter_V
    iter=iter+1;
    for k_ind = 1:N
        for z_ind = 1:M
            for kp_ind = 1:N
                Vtemp(k_ind,z_ind,kp_ind) = u(k_ind,z_ind,kp_ind) + beta* (P(z_ind,:)*V(kp_ind,:)');
            end
             [ Vnew(k_ind,z_ind) ,kp_ind ] = max(Vtemp(k_ind,z_ind,:));
             c_pol_disc(k_ind,z_ind) = c(k_ind,z_ind, kp_ind);
             L_pol_disc(k_ind,z_ind) = L(k_ind,z_ind, kp_ind);
             kprime_pol_disc(k_ind,z_ind) = kgrid(kp_ind);
             kprime_ind(k_ind, z_ind) = kp_ind;
        end
    end
    dV = max(abs(V-Vnew),[],'all');
    V=Vnew;
%     disp('dV');
%     disp(dV);
end

V_disc=V;
temp(1)=toc;

[x,y] = ndgrid(kgrid, exp(Z));
cprime_vfi_disc =  griddata(x,y, c_pol_disc, kprime_pol_disc(:), y(:)  ) ;
cprime_vfi_disc  = reshape(cprime_vfi_disc, [N,M]);
lprime_vfi_disc = griddata(x,y, L_pol_disc, kprime_pol_disc(:),y(:));
lprime_vfi_disc = reshape(lprime_vfi_disc, [N,M]  );

EE_error_vfi_disc = c_pol_disc.^(-gamma) - beta .*   cprime_vfi_disc.^(-gamma).*...
                  (exp(Z) .* alpha .* kprime_pol_disc .^(alpha-1) .* lprime_vfi_disc.^(1-alpha) + 1-delta   )  * P'   ;




% simulate stochastic transition
[xgrid,ygrid] = ndgrid(kgrid,exp(Z));
kprime_pol_disc_linInt = griddedInterpolant(xgrid,ygrid,kprime_pol_disc,'linear');
c_pol_disc_linInt = griddedInterpolant(xgrid,ygrid,c_pol_disc,'linear');

k_disc_sim =zeros(N_sim,T);
L_disc_sim =zeros(N_sim,T);
c_disc_sim =zeros(N_sim,T);
y_disc_sim =zeros(N_sim,T);
inv_disc_sim =zeros(N_sim,T);


for sim=1:N_sim
    k_disc_sim(sim,1)=kbar;
    for t=2:T
        % here I use interpolation of the policy function given z(sim,t-1) at k_disc_sim(sim,t-1). since
        % kprime is on the grid, this is not actually necessary
        % I also use interpolation of the policy function for L_disc_sim(sim,t-1) given z(sim,t-1) at k_disc_sim(sim,t-1)
        k_disc_sim(sim,t)=kprime_pol_disc_linInt(k_disc_sim(sim,t-1),Z_sim(sim,t-1));
        c_disc_sim(sim,t-1)=c_pol_disc_linInt(k_disc_sim(sim,t-1),Z_sim(sim,t-1));
    end
    c_disc_sim(sim,T)=c_pol_disc_linInt(k_disc_sim(sim,T),Z_sim(sim,T));
    L_disc_sim(sim,:) = ( 1/theta .* c_disc_sim(sim,:).^(-gamma) .* (1-alpha).* Z_sim(sim,:) .* k_disc_sim(sim,:).^(alpha)).^(1/(phi+alpha));   
    y_disc_sim(sim,:) = Z_sim(sim,:).*k_disc_sim(sim,:).^alpha.* L_disc_sim(sim,:).^(1-alpha);
    inv_disc_sim(sim,:) = y_disc_sim(sim,:) - c_disc_sim(sim,:) ;
end

%% 
% ===========================================
%  Option II:   Value function iteration with interpolation
% ===========================================

dV=1; 
    uu = zeros(N,M);
    uu = uu-Inf;
    [x,y]=ndgrid(kgrid,exp(Z));
    cc= x.^alpha.*1.1*Lbar.^(1-alpha) - delta.*x;
    ll= (1/theta .* cc.^(-gamma) .* (1-alpha) .* x.^(alpha)).^(1/(phi+alpha));
    uu= (cc.^(1-gamma)-1)./(1-gamma) - theta./(1+phi) .* ll.^(1+phi);
    V = uu./(1-beta);
    [x, y] = ndgrid(kgrid,exp(Z));
    V = griddedInterpolant(x,y,V);
    V_disc = uu./(1-beta);
% V = zeros(N,M);
% [x,y] = ndgrid(kgrid, exp(Z));
% V = griddedInterpolant(x,y,V);
% use previous solution as initial guess
iter=0;
tic
% a better guess: V=V_disc;
% initiate policies
c_interp=c_pol_disc; L_interp=L_pol_disc;
Vnew = zeros(N,M);
kprime_interp = zeros(N,M);
% loop over policy function
while dV>criter_V
    iter=iter+1;
    for k_ind=1:N
        for z_ind=1:M
            fun = @(kprime) -util(kprime,V ,gamma,theta, delta ,phi,beta, alpha, P, Z, kgrid ,k_ind, z_ind);
            [xx,yy]=fminsearch(fun, [cbar,Lbar]) ;
            Vnew(k_ind,z_ind)=util(xx,V ,gamma,theta, delta ,phi,beta, alpha, P, Z, kgrid,k_ind,z_ind)  ;
            c_interp(k_ind,z_ind)=xx(1);
            L_interp(k_ind,z_ind)=xx(2);
            kprime_interp(k_ind,z_ind)=exp(Z(z_ind))*kgrid(k_ind)^alpha*xx(2)^(1-alpha) + (1-delta)*kgrid(k_ind)-xx(1) ;
        end
    end
    if isreal(c_interp)==0 || isreal(L_interp)==-3
        keyboard
    end
    dV=max(max(abs(Vnew-V_disc)));
    [x, y] = ndgrid(kgrid,exp(Z));
    V = griddedInterpolant(x,y,Vnew,'linear');
    V_disc = V(x,y);
%     disp('dV')
%     disp(dV)
end
temp(2)=toc;

cprime_vfiIn =  griddata(x,y, c_interp, kprime_interp(:), y(:)  ) ;
cprime_vfiIn  = reshape(cprime_vfiIn, [N,M]);
lprime_vfiIn = griddata(x,y, L_interp, kprime_interp(:),y(:));
lprime_vfiIn = reshape(lprime_vfiIn, [N,M]  );

EE_error_vfiIn = c_interp.^(-gamma) - beta .*   cprime_vfiIn.^(-gamma).*...
                  (exp(Z) .* alpha .* kprime_interp .^(alpha-1) .* lprime_vfiIn.^(1-alpha) + 1-delta   )  * P'   ;

              
% simulate stochastic transition
[xgrid,ygrid] = ndgrid(kgrid,exp(Z));
kprime_pol_vfiIn_linInt = griddedInterpolant(xgrid,ygrid,kprime_interp,'linear');
c_pol_vfiIn_linInt = griddedInterpolant(xgrid,ygrid,c_interp,'linear');

k_vfiIn_sim =zeros(N_sim,T);
L_vfiIn_sim =zeros(N_sim,T);
c_vfiIn_sim =zeros(N_sim,T);
y_vfiIn_sim =zeros(N_sim,T);
inv_vfiIn_sim =zeros(N_sim,T);

for sim=1:N_sim
    k_vfiIn_sim(sim,1)=kbar;
    for t=2:T
        % here I use interpolation of the policy function given z(sim,t-1) at k_disc_sim(sim,t-1). since
        % kprime is on the grid, this is not actually necessary
        % I also use interpolation of the policy function for L_disc_sim(sim,t-1) given z(sim,t-1) at k_disc_sim(sim,t-1)
        k_vfiIn_sim(sim,t)=kprime_pol_vfiIn_linInt(k_vfiIn_sim(sim,t-1),Z_sim(sim,t-1));
        c_vfiIn_sim(sim,t-1)=c_pol_vfiIn_linInt(k_vfiIn_sim(sim,t-1),Z_sim(sim,t-1));
    end
    c_vfiIn_sim(sim,T)=c_pol_vfiIn_linInt(k_vfiIn_sim(sim,T),Z_sim(sim,T));
    L_vfiIn_sim(sim,:) = ( 1/theta .* c_vfiIn_sim(sim,:).^(-gamma) .* (1-alpha).* Z_sim(sim,:) .* k_vfiIn_sim(sim,:).^(alpha)).^(1/(phi+alpha));   
    y_vfiIn_sim(sim,:) = Z_sim(sim,:).*k_vfiIn_sim(sim,:).^alpha.* L_vfiIn_sim(sim,:).^(1-alpha);
    inv_vfiIn_sim(sim,:) = y_vfiIn_sim(sim,:) - c_vfiIn_sim(sim,:) ;
end





%% EGM
tic
[x,y] = ndgrid(kgrid, exp(Z));
c_policy_egm = y.*Lbar.^(1-alpha).*x.^alpha+(1-delta).*x - x; %initial guess
l_policy_egm  = (1/theta*c_policy_egm.^(-gamma).*(1-alpha).*y.*x.^alpha).^(1/(phi+alpha)); %l'(k',A')
k_egm = zeros(N,M); %k(k',A)
c_egm = zeros(N,M); %c(k',A)
l_egm = zeros(N,M);

% uprime = u^(-gamma)
error = 1;
iter = 0;
while error >=criter_V
    iter = iter  +1;
    for kp_ind =1:N
        for z_ind = 1:M
            c_egm(kp_ind,z_ind) = (  beta *  P(z_ind,:)  * ...
                                        (   c_policy_egm(kp_ind,:).^(-gamma).*...
                                          ( alpha.*exp(Z).*kgrid(kp_ind).^(alpha-1) .* l_policy_egm(kp_ind,:).^(1-alpha) ...
                                          + (1-delta) )  )'  )^(-1/gamma);
            % Find capital and labor in the current period.
           [xx,yy] = fsolve( @(x) FOC2(x,c_egm(kp_ind,z_ind),kgrid(kp_ind),exp(Z(z_ind)),gamma,theta,phi,alpha,beta,delta),[kbar,Lbar]  );
           k_egm(kp_ind,z_ind) = xx(1);
           l_egm(kp_ind,z_ind) = xx(2);
        end
    end
    temp = ones(N,1) * exp(Z);
    [x,y] = ndgrid(kgrid,exp(Z));
    new_c_policy_egm = griddata(k_egm(:),temp(:),c_egm(:),x,y,'v4') ;
    error = max(new_c_policy_egm-c_policy_egm,[],'all');
    
    c_policy_egm = new_c_policy_egm;
    l_policy_egm = ( 1/theta .* c_policy_egm.^(-gamma).*(1-alpha).*y.* x .^(alpha) ).^(1/(phi+alpha));
end

kprime_policy_egm = y.* x.^(alpha) .* l_policy_egm.^(1-alpha) + (1-delta) .* x - c_policy_egm;

% Euler Equation Error
kprime_policy_egm_int = griddedInterpolant(x,y, kprime_policy_egm);
l_policy_egm_int = griddedInterpolant(x,y, l_policy_egm);
c_policy_egm_int = griddedInterpolant(x,y, c_policy_egm);

EE_error_egm = zeros(N,M);
kprime_egm = kprime_policy_egm_int(x,y);
lprime_egm = l_policy_egm_int(kprime_egm,y);
cprime_egm = c_policy_egm_int(kprime_egm,y);
EE_error_egm = c_policy_egm.^(-gamma) - beta .*   cprime_egm.^(-gamma).*...
                  (exp(Z) .* alpha .* kprime_egm .^(alpha-1) .* lprime_egm.^(1-alpha) + 1-delta   )  * P'   ;

% simulate stochastic transition
[xgrid,ygrid] = ndgrid(kgrid,exp(Z));
kprime_pol_egm_linInt = griddedInterpolant(xgrid,ygrid,kprime_pol_disc,'linear');
c_pol_egm_linInt = griddedInterpolant(xgrid,ygrid,c_pol_disc,'linear');

k_egm_sim =zeros(N_sim,T);
L_egm_sim =zeros(N_sim,T);
c_egm_sim =zeros(N_sim,T);
y_egm_sim =zeros(N_sim,T);
inv_egm_sim =zeros(N_sim,T);

for sim=1:N_sim
    k_egm_sim(sim,1)=kbar;
    for t=2:T
        % here I use interpolation of the policy function given z(sim,t-1) at k_disc_sim(sim,t-1). since
        % kprime is on the grid, this is not actually necessary
        % I also use interpolation of the policy function for L_disc_sim(sim,t-1) given z(sim,t-1) at k_disc_sim(sim,t-1)
        k_egm_sim(sim,t)=kprime_pol_vfiIn_linInt(k_egm_sim(sim,t-1),Z_sim(sim,t-1));
        c_egm_sim(sim,t-1)=c_pol_vfiIn_linInt(k_egm_sim(sim,t-1),Z_sim(sim,t-1));
    end
    c_egm_sim(sim,T)=c_pol_vfiIn_linInt(k_egm_sim(sim,T),Z_sim(sim,T));
    L_egm_sim(sim,:) = ( 1/theta .* c_egm_sim(sim,:).^(-gamma) .* (1-alpha).* Z_sim(sim,:) .* k_egm_sim(sim,:).^(alpha)).^(1/(phi+alpha));   
    y_egm_sim(sim,:) = Z_sim(sim,:).*k_egm_sim(sim,:).^alpha.* L_egm_sim(sim,:).^(1-alpha);
    inv_egm_sim(sim,:) = y_egm_sim(sim,:) - c_egm_sim(sim,:) ;
end
toc





%%

% ==============
% Figures
% ==============

plotz=3;
figure('Name','Policy functions') % this is not asked
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 18]);
subplot(3,1,1)
title('Consumption','Interpreter','Latex','fontsize',13)
hold on
plot([c_pol_disc(:,plotz)])
plot([c_interp(:,plotz)])
plot([c_policy_egm(:,plotz)])
xlabel('Capital','Interpreter','Latex','fontsize',13)
subplot(3,1,2)
title('Labor supply','Interpreter','Latex','fontsize',13)
hold on
plot([L_pol_disc(:,plotz)])
plot([L_interp(:,plotz)])
plot([l_policy_egm(:,plotz)])
subplot(3,1,3)
title('K prime','Interpreter','Latex','fontsize',13)
hold on
plot([kprime_pol_disc(:,plotz)])
plot([kprime_interp(:,plotz)])
plot([kprime_egm(:,plotz)])
xlabel('Capital','Interpreter','Latex','fontsize',13)
h = legend('Discrete','Interpolate','EGM','Location', 'best','Orientation','Vertical');
set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'

saveas(gcf,".\figures\policy_functions.png")

% 3D policy function
figure('Name','Policy Function in 3D graph, VFI discrete')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 14]);
[x,y] = ndgrid(kgrid,exp(Z));
surf(x, y, kprime_pol_disc);
xlabel('$K_t$','Interpreter','Latex');
ylabel('$A_t$','Interpreter','Latex');
zlabel('$K_{t+1}$','Interpreter','Latex');
title('Policy Function, VFI Discrete');
saveas(gcf,".\figures\policy_func_vfi_disc.png")

% 3D policy function, VFI linear interpolation
figure('Name','Policy Function in 3D graph, VFI Interpolation')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 14]);
[x,y] = ndgrid(kgrid,exp(Z));
surf(x, y, kprime_interp);
xlabel('$K_t$','Interpreter','Latex');
ylabel('$A_t$','Interpreter','Latex');
zlabel('$K_{t+1}$','Interpreter','Latex');
title('Policy Function, VFI Discrete');
saveas(gcf,'.\figures\policy_func_vfi_int.png');

% 3D policy function, EGM
figure('Name','Policy Function in 3D graph, EGM')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 14]);
[x,y] = ndgrid(kgrid,exp(Z));
surf(x, y, kprime_egm);
xlabel('$K_t$','Interpreter','Latex');
ylabel('$A_t$','Interpreter','Latex');
zlabel('$K_{t+1}$','Interpreter','Latex');
title('Policy Function, VFI Discrete');
saveas(gcf,'.\figures\policy_func_egm.png');


% Euler Equation error
% 3D policy function
figure('Name','Euler Equation Errors in 3D graph, VFI discrete')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 14]);
[x,y] = ndgrid(kgrid,exp(Z));
surf(x, y, EE_error_vfi_disc*100);
xlabel('$K_t$','Interpreter','Latex');
ylabel('$A_t$','Interpreter','Latex');
zlabel('$Error$','Interpreter','Latex');
title('Euler Equation Errors(%), VFI Discrete');
saveas(gcf,".\figures\EEerror_vfi_disc.png")


% Euler Equation error
% 3D policy function
figure('Name','Euler Equation Errors in 3D graph, VFI Interpolation')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 14]);
[x,y] = ndgrid(kgrid,exp(Z));
surf(x, y, EE_error_vfiIn*100);
xlabel('$K_t$','Interpreter','Latex');
ylabel('$A_t$','Interpreter','Latex');
zlabel('$Errors$','Interpreter','Latex');
title('Euler Equation Errors(%), VFI Interpolation');
saveas(gcf,".\figures\EEerror_vfi_interp.png")




% Euler Equation error
% 3D policy function
figure('Name','Euler Equation Errors in 3D graph, VFI Interpolation')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 14]);
[x,y] = ndgrid(kgrid,exp(Z));
surf(x, y, EE_error_egm*100);
xlabel('$K_t$','Interpreter','Latex');
ylabel('$A_t$','Interpreter','Latex');
zlabel('Errors','Interpreter','Latex');
title('Euler Equation Errors(%), EGM');
saveas(gcf,".\figures\EEerror_egm.png")




figure('Name','Simulated time series')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 30]);
subplot(3,1,1)
hold on
title('Simulation of the capital stock')
%for sim=1:N_sim
sim=2
plot(k_disc_sim(sim,1:T),'k','Linewidth',1)
plot(k_vfiIn_sim(sim,1:T),'k:','Linewidth',1)
plot(k_egm_sim(sim,1:T),'k--','Linewidth',1)
%end
h = legend('DP discrete', 'DP interpol', 'DP EGM','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'
subplot(3,1,2)
hold on
title('Simulation of consumption')
%for sim=1:N_sim
sim=1;
plot(c_disc_sim(sim,1:T),'k','Linewidth',1)
plot(c_vfiIn_sim(sim,1:T),'k:','Linewidth',1)
plot(c_egm_sim(sim,1:T),'k--','Linewidth',1)
h = legend('DP discrete', 'DP interpol','DP EDM','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'
subplot(3,1,3)
hold on
title('Simulation of Labor')
%for sim=1:N_sim
sim=1
plot(L_disc_sim(sim,1:T),'k','Linewidth',1)
plot(L_vfiIn_sim(sim,1:T),'k:','Linewidth',1)
plot(L_egm_sim(sim,1:T),'k--','Linewidth',1)
h = legend('DP discrete', 'DP interpol', 'DP EDM', 'Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'
saveas(gcf,".\figures\simulation.png")




for i=1:N_sim
    rho_emp(i)=corr(log(Z_sim(i,1:end-1))',log(Z_sim(i,2:end))');
    rho_emp_cont(i)=corr(z_cont(i,1:end-1)',z_cont(i,2:end)');
    corr_Cy_disc(i)=corr(y_disc_sim(i,:)',c_disc_sim(i,:)');
    corr_invy_disc(i)=corr(y_disc_sim(i,:)',inv_disc_sim(i,:)');
    corr_Ly_disc(i)=corr(y_disc_sim(i,:)',L_disc_sim(i,:)');
    corr_Cy_interp(i)=corr(y_vfiIn_sim(i,:)',c_vfiIn_sim(i,:)');
    corr_invy_interp(i)=corr(y_vfiIn_sim(i,:)',inv_vfiIn_sim(i,:)');
    corr_Ly_interp(i)=corr(y_vfiIn_sim(i,:)',L_vfiIn_sim(i,:)');
        corr_Cy_egm(i)=corr(y_egm_sim(i,:)',c_egm_sim(i,:)');
    corr_invy_egm(i)=corr(y_egm_sim(i,:)',inv_egm_sim(i,:)');
    corr_Ly_egm(i)=corr(y_egm_sim(i,:)',L_egm_sim(i,:)');
end
disp('Autocorrelation of Z discrete, continuous, theoretical')
disp([mean(rho_emp), mean(rho_emp_cont), rho    ])
disp('Standard devations of Z discrete, continuous, theoretical')
disp([mean(std(log(Z_sim)')),mean(std(z_cont')),(sigmaepsilon^2/(1-rho^2))^0.5])

disp('Standard Deviation of Y, C,I,L')
disp('VFI discrete')
disp([mean(std(y_disc_sim')), mean(std(c_disc_sim')), mean(std(inv_disc_sim')), mean(std(L_disc_sim'))])
disp('VFI interpolation')
disp([mean(std(y_vfiIn_sim')), mean(std(c_vfiIn_sim')), mean(std(inv_vfiIn_sim')), mean(std(L_vfiIn_sim'))])
disp('EGM')
disp([mean(std(y_egm_sim')), mean(std(c_egm_sim')), mean(std(inv_egm_sim')), mean(std(L_egm_sim'))])


disp('Correlation with y of C,I,L')
disp('VFI discrete')
disp([mean(corr_Cy_disc),mean(corr_invy_disc),mean(corr_Ly_disc)])
disp('VFI interpolation')
disp([mean(corr_Cy_interp),mean(corr_invy_interp),mean(corr_Ly_interp)])
disp('EGM')
disp([mean(corr_Cy_egm),mean(corr_invy_egm),mean(corr_Ly_egm)])

