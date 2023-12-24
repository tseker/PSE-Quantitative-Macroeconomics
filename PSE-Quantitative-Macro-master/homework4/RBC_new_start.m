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

% ============
% options and convergence criteria
% ============
criter_V = 1e-6; % conv criterion for value function
N=80; % number of grid points
linear=1; % grid linear or not
M=3; % number of support points for the shock
T=1000; % periods for transition
N_sim=100; % number of simulations

% ==============
% Problem I: Calibration
% ==============

beta=1/Rbar;
delta=InvKrat;
alpha=Kshare;
margprod=1/beta-1+delta;
KNrat=(margprod/alpha)^(1/(alpha-1));
wbar=(1-alpha)*(margprod/alpha)^(alpha/(alpha-1));
cKrat = margprod/alpha-delta;
theta = wbar/(   Lbar^(gamma+psi)* ( (margprod/alpha) -delta)^gamma * (margprod/alpha)^(gamma/(alpha-1)) );
kbar = KNrat*Lbar;

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
%z= zeros(N,1);
rng(1);
z = zeros(N_sim, T);
for j=1:N_sim
    z(j,1)=3;
    for t=2:T
        z(j,t)=sum( cumsum(P(z(j,t-1),:))<rand(1,M)  )+1;
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
Z_sim=Z_lev(z);
% cbar is not the same




%%


% =====================
% Option III:   Log-linearization, simulation of 1000 periods
% =====================


% some ratios as function of parameters
ybar=kbar^alpha*Lbar^(1-alpha);
cbar=ybar-delta*kbar;
ckrat=cbar/kbar;
Rbar=alpha*(ybar/kbar);

% a. write system as A_1 E[x_t+1]=A_2  x_t+B_1 epsilon_t
% order k,a,c,l
a= [beta*Rbar*(1-alpha), -beta*Rbar, gamma, -beta*Rbar*(1-alpha); 0 0 0 0; kbar/ybar, 0, 0, 0; 0 1 0 0];

b=	[0 0 gamma 0; alpha 1 -gamma -(alpha+psi); 1/beta*kbar/ybar 1 -cbar/ybar (1-alpha); 0 rho 0 0];

% re-order to have state variables first
% order k,a,c,n
nk=2;
[f,p] = solab(a,b,nk);

% extract cons and lab policies
c_polfunc=f(1,:);
n_polfunc=f(2,:);
LOM=p(1,:);

% policy functions are in log deviations. so need to pre-multiply by cbar
% and add cbar to get levels
c_pol_lin = zeros(length(kgrid),M);
n_pol_lin = zeros(length(kgrid),M);
k_pol_lin = zeros(length(kgrid),M);
for j=1:M
    c_pol_lin(:,j)=cbar*exp(c_polfunc(1)*log(kgrid/kbar)+c_polfunc(2)*Z(j));
    n_pol_lin(:,j)=Lbar*exp(n_polfunc(1)*log(kgrid/kbar)+n_polfunc(2)*Z(j));
    k_pol_lin(:,j)=kbar*exp(LOM(1)*log(kgrid/kbar)+LOM(2)*Z(j)); 
end


k_sim_lin=zeros(N_sim,T);
k_sim_lin(:,1) = kbar*ones(N_sim,1);
c_sim_lin=zeros(N_sim,T);
L_sim_lin = zeros(N_sim,T);
inv_sim_lin = zeros(N_sim,T);
y_sim_lin = zeros(N_sim, T);


% generate a 

for j=1:N_sim
    for t=2:T
        c_sim_lin(j,t-1) = cbar*exp(c_polfunc(1)*log(k_sim_lin(j,t-1)/kbar)+c_polfunc(2)*z_cont(j,t-1));
        L_sim_lin(j,t-1) = Lbar*exp(n_polfunc(1)*log(k_sim_lin(j,t-1)/kbar)+n_polfunc(2)*z_cont(j,t-1));
        k_sim_lin(j,t) = kbar*exp(LOM(1)*log(k_sim_lin(j,t-1)/kbar)+LOM(2)*z_cont(j,t-1));
        y_sim_lin(j,t-1) = exp(z_cont(j,t-1))*(k_sim_lin(j,t-1)^alpha)*L_sim_lin(j,t-1)^(1-alpha);
        inv_sim_lin(j,t-1) = k_sim_lin(j,t)-(1-delta)*k_sim_lin(j,t-1);
    end
    
    c_sim_lin(j,T)=cbar;
    L_sim_lin(j,T)=Lbar;
    inv_sim_lin(j,T)=delta*kbar;
    y_sim_lin(j,T)=exp(z_cont(j,T-1))*(kbar^alpha)*(Lbar^(1-alpha));
end

k_sim_lin_discrete=zeros(N_sim,T);
k_sim_lin_discrete(:,1) = kbar*ones(N_sim,1);
c_sim_lin_discrete=zeros(N_sim,T);
L_sim_lin_discrete = zeros(N_sim,T);
inv_sim_lin_discrete = zeros(N_sim,T);
y_sim_lin_discrete = zeros(N_sim, T);



for j=1:N_sim
    for t=2:T
        c_sim_lin_discrete(j,t-1) = cbar*exp(c_polfunc(1)*log(k_sim_lin_discrete(j,t-1)/kbar)+c_polfunc(2)*log(Z_sim(j,t-1)));
        L_sim_lin_discrete(j,t-1) = Lbar*exp(n_polfunc(1)*log(k_sim_lin_discrete(j,t-1)/kbar)+n_polfunc(2)*log(Z_sim(j,t-1)));
        k_sim_lin_discrete(j,t) = kbar*exp(LOM(1)*log(k_sim_lin_discrete(j,t-1)/kbar)+LOM(2)*log(Z_sim(j,t-1))   );
        y_sim_lin_discrete(j,t-1) = Z_sim(j,t-1)*(k_sim_lin_discrete(j,t-1)^alpha)*L_sim_lin(j,t-1)^(1-alpha);
        inv_sim_lin_discrete(j,t-1) = y_sim_lin_discrete(j,t-1)-c_sim_lin_discrete(j,t-1);
    end
    
    c_sim_lin_discrete(j,T)=cbar;
    L_sim_lin_discrete(j,T)=Lbar;
    inv_sim_lin_discrete(j,T)=delta*kbar;
    y_sim_lin_discrete(j,T)=Z_sim(j,T-1)*(kbar^alpha)*(Lbar^(1-alpha));
end


%%

% ==============
% Option IV  Transition after one-time shock
% ==============

T = 150;

% initialise vectors
options=optimset('MaxIter',1000,'MaxFunEval',10000);
mu_trans = zeros(T,1);
mu_trans(1)=sigmaepsilon;
for t=2:T
    mu_trans(t)=rho*mu_trans(t-1);
end
z_trans=exp(mu_trans); clear mu_trans

k_0=kbar;
x0=[kbar;cbar;Lbar]*ones(1,T).*1.1; x0=x0';
x0(1,1)=k_0;
f=@(x) rbc_obj_endog_N(x,alpha,beta,gamma,delta,theta,psi,k_0,cbar,kbar,Lbar,z_trans);
[x,fsol]=fsolve(f,x0,options);

k_trans=x(:,1);
c_trans=x(:,2);
n_trans=x(:,3);
y_trans= z_trans.*k_trans.^alpha.* n_trans.^(1-alpha);
i_trans= y_trans-c_trans;
r_trans=  z_trans.*alpha.*k_trans.^(alpha-1).* n_trans.^(1-alpha);
w_trans=  z_trans.*(1-alpha).*k_trans.^alpha.* n_trans.^(-alpha);

T = 1000;

% for j=1:N_sim
%     for t=1:T
%         if t>1
%             k_det_sim(j,t:t+T-1)=k_det_sim(j,t:t+T-1)+(k_trans'-kbar)*(z_cont(j,t)/z_cont(j,t-1)^rho-1)/(sigmaepsilon);
%             c_det_sim(j,t:t+T-1)=[Fill this in];
%             L_det_sim(j,t:t+T-1)=[Fill this in];
%             y_det_sim(j,t-1)=[Fill this in];
%         else
%             k_det_sim(j,t:t+T-1)=[Fill this in];
%             c_det_sim(j,t:t+T-1)=[Fill this in];
%             L_det_sim(j,t:t+T-1)=[Fill this in];
%             y_det_sim(j,T)=[Fill this in];
%         end
%         %k_det_sim(j,t+1:t+T-1)=k_det_sim(j,t+1:t+T-1)+kprime_trans(1:end-1)'*Z_sim(j,t);
%     end
%     
%     for t=2:T
%         inv_det_sim(j,t-1)=k_det_sim(j,t)-(1-delta)*k_det_sim(j,t-1);
%     end
%     inv_det_sim(j,T)=Z_sim(j,T)*k_det_sim(j,T)^alpha*L_det_sim(j,T)^(1-alpha)-c_det_sim(j,T);
% end
% c_det_sim=c_det_sim(:,1:T);
% L_det_sim=L_det_sim(:,1:T);
% k_det_sim=k_det_sim(:,1:T);

%%
plotting


