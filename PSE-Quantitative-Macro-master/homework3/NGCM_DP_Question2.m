% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all
% change location to the folder where this m file is saved
mfile_name          = mfilename('fullpath');  % This line of code retrieves the full path of the current MATLAB file and assigns it to the variable mfile_name. The mfilename function with the argument 'fullpath' returns the full path of the currently executing file, including the file name and extension. The full path is then assigned to the variable mfile_name.
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% ============
% parameters
% ============
theta = 0.4; % capital share
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
gamma = 2.0000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta = 0.1;

% ============
% options and convergence criteria
% ============

Howard =0; % set to 1 if you want to do policy fct iteration / Howard improvement
criter_V = 1e-7; % conv criterion for value function
N=50; % number of grid points
linear=1; % grid linear or not
T=250;
%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(theta))^(1/(theta-1));

% ==============
% 0. Grids,  etc
% ==============

% center the grid around mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(theta))^(1/(theta-1))
% the grid
if delta==1
    if linear==1
        kgrid=linspace(kbar/2,2*kbar,N);
    else
        temp=linspace(0,0.5,N).^5/0.5^5*(2*kbar-kbar/2);
        kgrid=kbar/2+temp
    end
else
    kgrid=linspace(kbar*3/4,kbar*5/4,N);
end

% ==============
% 1. analytical case delta=1, finite T and infinite T
% ==============

% a. analytical policies, finite and infinite horizon

% ==============
% 2. Value function iteration (solve with and w/o Howard / policy function iteration), infinite T
% ==============

% one period return
for i=1:N
    for j=1:N
        c(i,j)= kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j);
        if c(i,j)>0
            u(i,j)= (c(i,j)^(1-gamma)-1)/(1-gamma);
        else
            u(i,j)=-1e50*((kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j))<=0);
        end
    end
end

% this give a consumption and utility vector when full depreciation
for i=1:N
    cc(i)= kgrid(i)^theta-delta*kgrid(i);
    if cc(i)>0
        uu(i)= (cc(i)^(1-gamma)-1)/(1-gamma);
    else
        uu(i)=-1e50*((kgrid(i)^theta-delta*kgrid(i))<=0);
    end
end
% initial guesses
dV=1;
V= zeros(N,1);% 0 works but not efficient 

%stationary linear consumption rule: a good guess is
%u(k^theta-delta*k)/(1-beta)
% for i=1:N
% V(i)=uu(i)/(1-beta)
VB= zeros(N,1);
for i=1:N
    VB(i)= uu(i)/(1-beta);
end
kprime = zeros(N,1);
max_index= zeros(N,1); 
iter=0;
tic
V=VB;

% main loop
while dV>criter_V
    iter=iter+1;
    for i=1:N % loop over capital today
        for j=1:N %... and over capital tomorrow
            VV(i,j)= u(i,j)+ beta*V(j); % calculate value for each i,j
        end
        % take maximum over capital tomorrow
        [Vnew(i),max_index(i)] = max(VV(i,:));
        % record policy function
        kprime(i) = kgrid(max_index(i));
    end
    % Howard - doesn't help much here
    if Howard==1 && iter>0
        dVV=1;
        while dVV>criter_V
            for i=1:N
                Vnewhow(i)=u(i,max_index(i))+beta*Vnew(max_index(i));
                clear temp
            end
            dVV=max(max(abs(Vnewhow-Vnew)));
            %disp(dVV)
            Vnew=Vnewhow;
        end
    end   
    % calculate convergence criterion
    dV=max(max(abs(Vnew-V)));
    % updated value function
    V=Vnew;
    disp('dV')
    disp(dV)
end

V_disc_VFI=V;
toc
 
kcompare=zeros(N,1);
for i=1:N
    kcompare(i)= beta*theta*(kgrid(i)^theta);
end
vcompare=zeros(N,1);

% Create a figure
figure;
% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(kgrid,kprime,'Linewidth',0.9);
hold on;
plot(kgrid,kcompare,'Linewidth',0.9);
title('Deterministic VFI Policy Function $k_{t+1}(k_{t})$, Nonlinear, $N_k=200$, $\delta=1$, $\sigma=1$','interpreter','latex','FontSize',12);
xlabel('Capital $k_t$','interpreter','latex','FontSize',16);
ylabel('Capital $k_{t+1}$','interpreter','latex','FontSize',16);
xline(kbar);
yline(kbar);
% Plot the 45-degree line, help to check the steady state point
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k^*'; % The label for the special point
text(0.5, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(kgrid,V_disc_VFI,'Linewidth',0.9);
title('Deterministic VFI Value Function $V_{t}(k_{t})$, Nonlinear, $N_k=200$, $\delta=1$, $\sigma=1$','interpreter','latex','FontSize',12);
xlabel('Capital $k_t$','interpreter','latex','FontSize',16);
ylabel('Maximized Utility $V_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'.\figures\rc_sigmadelta1_nonlinear.png');
hold off; 

% Euler Equation Error
% consumption vector today
c1=(kgrid.^theta+(1-delta)*kgrid)'-kprime;
%consumption vector at choice kprime tomorrow
c2=kprime.^theta+(1-delta)*kprime-kprime(max_index);
%marginal productivity
margprod=theta*kprime.^(theta-1)+1-delta;
%implied c1
c1_imp=(beta*margprod).^(-1/gamma).*c2;
%Error
EEerror_disc=(c1-c1_imp)./c1_imp*100;%percentage deviation 
maxEEerror_disc=max(abs(EEerror_disc));
%plot the figure
figure(2);
hold on; 
plot(kgrid,EEerror_disc,'Linewidth',0.9);
title('Euler Equation Errors, Baseline, $N_k=50$','interpreter','latex','FontSize',12);
xlabel('Capital $k_t$','interpreter','latex','FontSize',16);
ylabel('Euler Equation Error (\%)','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'.\figures\eee_linear_50.png');
hold off;