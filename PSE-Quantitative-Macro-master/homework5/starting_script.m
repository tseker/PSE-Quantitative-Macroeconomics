
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
psi=1; %inverse Frisch elasticity


% Calibration
beta=1/Rbar;
delta=InvKrat;
alpha=Kshare;
margprod=1/beta-1+delta;
KNrat=(margprod/alpha)^(1/(alpha-1));
wbar=(1-alpha)*(margprod/alpha)^(alpha/(alpha-1));
cKrat = margprod/alpha-delta;
theta = wbar/(   Lbar^(gamma+psi)* ( (margprod/alpha) -delta)^gamma * (margprod/alpha)^(gamma/(alpha-1)) );
kbar = KNrat*Lbar;



% ============
% options and convergence criteria
% ============
criter_V = 1e-6; % conv criterion for value function
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
for j=1:N_sim
    z(j,1)=3;
    for t=2:T
        z(j,t)=[Fill this in];
    end
end

% simulate continuous Markov Process
for j=1:N_sim
    z_cont(j,1)=1;
    for t=2:T
        z_cont(j,t)=[Fill this in ];
    end
end

Z_lev=exp(Z);
Z_sim=Z_lev(z);



% ==============
% Option I:  Fully discretized value function iteration
% ==============
options=optimset('MaxIter',1000,'MaxFunEval',1000);
toc
% one period return
for i=1:N
    for k=1:M
        for j=1:N
            if kgrid(j)>exp(Z(k))*kgrid(i)^alpha*1^(1-alpha)+(1-delta)*kgrid(i)
                c(i,j,k)=-1;
                L(i,j,k)=1;
            else
                % here you need to solve the intratemporal FOC for L
                [xx,fval,exitflag]=fsolve(@(x) FOC(x,kgrid(i),Z(k),gamma,theta,psi,alpha,beta,delta,kgrid(j)),[Lbar],options);
                exit(i,j,k)=exitflag;
                L(i,j,k)=xx;
                c(i,j,k)=[Fill this in];
                
            end
            if c(i,j,k)<=0 || exit(i,j,k)~=1
                u(i,j,k)=-1e10;
            else
                if L(i,j,k)>=1
                    L(i,j,k)=1;
                    c(i,j,k)=exp(Z(k))*kgrid(i)^alpha*L(i,j,k)^(1-alpha)+(1-delta)*kgrid(i)-kgrid(j);
                elseif  L(i,j,k)<=0
                    L(i,j,k)=0;
                    c(i,j,k)=exp(Z(k))*kgrid(i)^alpha*L(i,j,k)^(1-alpha)+(1-delta)*kgrid(i)-kgrid(j);
                end
                u(i,j,k)=[Fill this in];
            end
        end
    end
end
% iterate on value function
dV=1;
Ytemp=(exp(Z)*kgrid).^alpha*(0.5)^(1-alpha);
% initial guess
V=1/(1-beta)*(((0.3*(Ytemp+(1-delta).*(ones(M,1)*kgrid))).^(1-gamma)-1)./(1-gamma)-theta*(0.3)^(1+psi)/(1+psi));
V=V';
iter=0;
tic

% iterate on Value Function
while dV>criter_V
    iter=iter+1;
    for i=1:N
        for k=1:M
            for j=1:N
                Vtemp(i,j,k)= u(i,j,k)+beta*P(k,:)*V(j,:)';
            end
            [Vnew(i,k),indic(i,k)]=[Fill this in];
            c_disc(i,k)=c(i,indic(i,k),k);
            L_disc(i,k)=L(i,indic(i,k),k);
            u_disc(i,k)=u(i,indic(i,k),k);
            kprime_disc(i,k)=kgrid(indic(i,k));
        end
    end
    dV=max(max(abs(Vnew-V)));
    V=Vnew;
    disp('dV')
    disp(dV)
end

V_disc=V;
temp(1)=toc;

% simulate stochastic transition
for sim=1:N_sim
    k_disc_sim(sim,1)=kbar;
    for t=2:T
        % here I use interpolation of the policy function given z(sim,t-1) at k_disc_sim(sim,t-1). since
        % kprime is on the grid, this is not actually necessary
        % I also use interpolation of the policy function for L_disc_sim(sim,t-1) given z(sim,t-1) at k_disc_sim(sim,t-1)
        k_disc_sim(sim,t)=interp1(kgrid,kprime_disc(:,z(sim,t-1)),k_disc_sim(sim,t-1));
        inv_disc_sim(sim,t-1)=[Fill this in];
        L_disc_sim(sim,t-1)=[Fill this in];
        c_disc_sim(sim,t-1)=[Fill this in];
        y_disc_sim(sim,t-1)=[Fill this in];

    end
    inv_disc_sim(sim,T)=delta*kbar;
    L_disc_sim(sim,T)=[Fill this in];
    c_disc_sim(sim,T)=Z_sim(sim,T)*k_disc_sim(sim,T)^alpha*L_disc_sim(sim,T)^(1-alpha)+(1-delta)*k_disc_sim(sim,T)-interp1(kgrid,kprime_disc(:,z(sim,T)),k_disc_sim(sim,T));
        y_disc_sim(sim,T)=Z_lev(z(sim,T))*k_disc_sim(sim,T)^alpha*L_disc_sim(sim,T)^(1-alpha);

end

% ===========================================
%  Option II:   Value function iteration with interpolation
% ===========================================
dV=1; clear cons clear L
Ytemp=(exp(Z)*kgrid).^alpha*(0.5)^(1-alpha);
% use previous solution as initial guess
V=1/(1-beta)*((c_disc).^(1-gamma)-1)./(1-gamma)-theta*L_disc.^(1+psi)/(1+psi);
iter=0;
tic
% a better guess: V=V_disc;
% initiate policies
c_interp=c_disc; L_interp=L_disc;
% loop over policy function
while dV>criter_V
    iter=iter+1;
    
    for i=1:N
        for k=1:M
            [xx,yy]=fminsearch(@(x) -util([Fill this in]),[c_interp(i,k);L_interp(i,k)],options);
            Vnew(i,k)=util([Fill this in]);
            c_interp(i,k)=xx(1);
            L_interp(i,k)=xx(2);
            kprime_interp(i,k)=exp(Z(k))*kgrid(i)^alpha*L_interp(i,k)^(1-alpha)+(1-delta)*kgrid(i)-c_interp(i,k);
            %Vnew(i,k)=-Vnew(i,k);
        end
    end
    if isreal(c_interp)==0 || isreal(L_interp)==0
        keyboard
    end
    dV=max(max(abs(Vnew-V)));
    V=Vnew;
    disp('dV')
    disp(dV)
end
temp(2)=toc;


% simulate stochastic transition
for sim=1:N_sim
    k_interp_sim(sim,1)=kbar;
    for t=2:T
        k_interp_sim(sim,t)=interp1(kgrid,kprime_interp(:,z(sim,t-1)),k_interp_sim(sim,t-1));
        inv_interp_sim(sim,t-1)=[Fill this in];
        c_interp_sim(sim,t-1)=[Fill this in];
        L_interp_sim(sim,t-1)=[Fill this in];
        y_interp_sim(sim,t-1)=[Fill this in];
    end
   
    xx=interp1(kgrid,kprime_interp(:,z(sim,t-1)),k_interp_sim(sim,t-1));
     inv_interp_sim(sim,T)=[Fill this in];
    c_interp_sim(sim,T)=[Fill this in];
    L_interp_sim(sim,T)=[Fill this in];
    y_interp_sim(sim,T)=[Fill this in];
end

% EGM
c_policy = Z.*lbar.^(1-alpha).*kgrid.^alpha+(1-delta)*kgrid - kgrid; %initial guess
k_t = zeros(M,N); RHS = k_t; M = k_t; c_new = k_t; c_policy_new = k_t; %initialization
err = 1;
while err > 1e-5
    for j = 1:M
        for i = 1:N
            RHS(j,i) = 0;
            for k = 1:M
                l = [xxx fill this in xxx]; %labor supply in t+1 for a given z_(t+1),k_(t+1)
                RHS(j,i) = [xxx fill this in xxx];
            end
            c_new(j,i) = [xxx fill this in xxx];
            M(j,i) = c_new(j,i)+kgrid(i);
            k_t(j,i) = newton(kgrid(i),c_new(j,i),Z(j),alpha,theta,mu,delta,sigma,M(j,i));
        end
        c_policy_new(j,:) = [xxx fill this in xxx];
    end
    err = sum(sum(abs(c_policy-c_policy_new)))
    c_policy = c_policy_new;
end
l_policy = [xxx fill this in xxx];
k_policy = [xxx fill this in xxx];


% ==============
% Figures
% ==============

plotz=3;
figure('Name','Policy functions') % this is not asked
subplot(3,1,1)
title('Consumption','Interpreter','Latex','fontsize',13)
hold on
plot([c_disc(:,plotz)])
plot([c_interp(:,plotz)])
plot([c_pol_lin(:,plotz)])
xlabel('Capital','Interpreter','Latex','fontsize',13)
subplot(3,1,2)
title('Labor supply','Interpreter','Latex','fontsize',13)
hold on
plot([L_disc(:,plotz)])
plot([L_interp(:,plotz)])
plot([n_pol_lin(:,plotz)])
subplot(3,1,3)
title('K prime','Interpreter','Latex','fontsize',13)
hold on
plot([kprime_disc(:,plotz)])
plot([kprime_interp(:,plotz)])
plot([k_pol_lin(:,plotz)])
xlabel('Capital','Interpreter','Latex','fontsize',13)
h = legend('Discrete','Interpolate','Linearised','Location', 'best','Orientation','Vertical');
set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'



figure('Name','Simulated time series')
subplot(3,1,1)
hold on
title('Simulation of the capital stock')
%for sim=1:N_sim
sim=1
plot(k_disc_sim(sim,1:T),'k','Linewidth',1)
plot(k_interp_sim(sim,1:T),'k:','Linewidth',1)
plot(k_sim_lin(sim,1:T),'k--','Linewidth',1)
plot(k_det_sim(sim,1:T),'k-+','Linewidth',1)
%end
h = legend('DP discrete', 'DP interpol', 'Linear','cumulative MIT shocks','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'
subplot(3,1,2)
hold on
title('Simulation of consumption')
%for sim=1:N_sim
sim=1
plot(c_disc_sim(sim,1:T),'k','Linewidth',1)
plot(c_interp_sim(sim,1:T),'k:','Linewidth',1)
plot(c_sim_lin(sim,1:T),'k--','Linewidth',1)
plot(c_det_sim(sim,1:T),'k-+','Linewidth',1)
h = legend('DP discrete', 'DP interpol', 'Linear','cumulative MIT shocks','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'
subplot(3,1,3)
hold on
title('Simulation of consumption')
%for sim=1:N_sim
sim=1
plot(L_disc_sim(sim,1:T),'k','Linewidth',1)
plot(L_interp_sim(sim,1:T),'k:','Linewidth',1)
plot(L_sim_lin(sim,1:T),'k--','Linewidth',1)
plot(L_det_sim(sim,1:T),'k-+','Linewidth',1)
h = legend('DP discrete', 'DP interpol', 'Linear','cumulative MIT shocks','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'



disp('Standard devations of Z discrete, continuous, theoretical')
disp([mean(std(log(Z_sim)')),mean(std(log(z_cont)')),(sigmaepsilon^2/(1-rho^2))^0.5])
for i=1:N_sim
    rho_emp(i)=corr(log(Z_sim(i,1:end-1))',log(Z_sim(i,2:end))');
    rho_emp_cont(i)=corr(log(z_cont(i,1:end-1))',log(z_cont(i,2:end))');
    corr_Cy_disc(i)=corr(y_disc_sim(i,:)',c_disc_sim(i,:)');
    corr_invy_disc(i)=corr(y_disc_sim(i,:)',inv_disc_sim(i,:)');
    corr_Ly_disc(i)=corr(y_disc_sim(i,:)',L_disc_sim(i,:)');
    corr_Cy_interp(i)=corr(y_interp_sim(i,:)',c_interp_sim(i,:)');
    corr_invy_interp(i)=corr(y_interp_sim(i,:)',inv_interp_sim(i,:)');
    corr_Ly_interp(i)=corr(y_interp_sim(i,:)',L_interp_sim(i,:)');
        corr_Cy_det(i)=corr(y_det_sim(i,:)',c_det_sim(i,:)');
    corr_invy_det(i)=corr(y_det_sim(i,:)',inv_det_sim(i,:)');
    corr_Ly_det(i)=corr(y_det_sim(i,:)',L_det_sim(i,:)');
        corr_Cy_lin(i)=corr(y_sim_lin(i,:)',c_sim_lin(i,:)');
    corr_invy_lin(i)=corr(y_sim_lin(i,:)',inv_sim_lin(i,:)');
    corr_Ly_lin(i)=corr(y_sim_lin(i,:)',L_sim_lin(i,:)');
end
disp('Autocorrelation of Z discrete, continuous')
disp([Fill this in])
disp('Standard devations of Z discrete, continuous')
disp([Fill this in])

disp('Correlation with y of C,I,L')
disp('VFI discrete')
disp([Fill this in])
disp('VFI interpolation')
disp([Fill this in])
disp('Linear')
disp([Fill this in])
disp('IRF')
disp([Fill this in])
disp('Standard devations of L,i,c relative to y')
disp('VFI discrete')
disp([Fill this in])
disp('VFI interpolation')
disp([Fill this in])
disp('Linear')
disp([Fill this in])
disp('IRF')
disp([Fill this in])
