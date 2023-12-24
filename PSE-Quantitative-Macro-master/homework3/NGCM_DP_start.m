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
theta=0.4; % capital share
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 2; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.1;

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

%%
% ==============
% 0. Grids,  etc
% ==============

% center the grid around mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(theta))^(1/(theta-1));

%Initial level of capital in the transition
k_0=kbar*0.25;

% the grid
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

%%
% ==============
% 1. analytical case delta=1, sigma = 1, finite T and infinite T
% ==============

% a. analytical policies, finite and infinite horizon
% V = alpha0 + alphak log(k_{t+1})
% k_prime = beta* theta * k^theta 

alphak = theta/(1-theta*beta);
alpha0 = log(1-theta*beta)/(1-beta)+ 1/(1-beta)* beta* theta /(1-beta*theta)*log(beta*theta);
% kbar_log
kbar_log=((1/beta)/(theta))^(1/(theta-1));
kgrid_log=linspace(kbar_log/2,2*kbar_log,N); % kgrid for log utilty
k_policy_analytical = beta* theta * kgrid_log.^theta; % analytical solution
[c_policy_log_vfi, k_policy_log_vfi, Vlog ]=  ValueFunIteration_Discrete( N, theta, 1, beta, 1, criter_V,0 );
policyk_error= k_policy_analytical - k_policy_log_vfi;

figure_num = 1;
figure(figure_num)
plot(kgrid_log, k_policy_analytical);
hold on
plot(kgrid_log, k_policy_log_vfi);
hold off;


% consumption vector today
for i = 1:length(k_policy_log_vfi)
    element = k_policy_log_vfi(i);
    index = find(kgrid_log == element);  % 找到匹配元素的位置
    positions(i) = index;
end

c1=(kgrid_log.^theta )' - k_policy_log_vfi;

% consumption vector at choice kprime tomorrow
c2=k_policy_log_vfi.^theta -k_policy_log_vfi(positions);

% marginal productivity
margprod = theta*k_policy_log_vfi.^(theta-1);

% implied c1
c1_imp = (beta*margprod).^(-1) .* c2;

%error
EEerror_log=(c1 - c1_imp)./c1_imp * 100; %percentage deviation
maxEEerror_disc = max(abs(EEerror_log));


fsize=16;

%Plot for error
    figure_num = figure_num+1;
    figure(figure_num)
    hold on
        plot(kgrid_log,EEerror_log)
        grid on
        xlabel('Stock of capital')
        ylabel('% error in Euler Equation')
        set(gca,'FontSize',fsize)
    hold off






%%

% ==============
% 2. Value function iteration (solve with and w/o Howard / policy function iteration), infinite T
% ==============

% iii. one period return
c = zeros(N,N);
u = zeros(N,N);

for i=1:N
    for j=1:N
        c(i,j)= kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j); %Feasibility constraint
        if c(i,j)>0 
            u(i,j)=(c(i,j)^(1-sigma)-1)/(1-sigma); %Utility function
        else
            u(i,j)=-1e50*((kgrid(i)^theta+(1-delta)*kgrid(i)-kgrid(j))<=0);
        end
    end
end

% iv. initial guesses
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
            if kgrid(j) >= 0 && kgrid(j) <= (1-delta)*kgrid(i)+kgrid(i)^theta
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
%     disp('dV')
%     disp(dV)
end

display(iter)

V_disc_VFI=V;
toc

%c. Policy function plot
    % Options
    lwidth=1.5;
    fsize=16;
    % Graph
    figure_num = figure_num+1;
    figure(figure_num)
        hold on
            plot(kgrid, kprime_vf, 'LineWidth', lwidth)
            grid on
            xlabel('Current stock of capital')
            ylabel('Next period''s stock of capital')
            set(gca,'FontSize',fsize)
            %title('Policy functions')
        hold off

% Euler equation errors in percent

% consumption vector today
c1=(kgrid.^theta + (1-delta)*kgrid)' - kprime_vf;

% consumption vector at choice kprime tomorrow
c2=kprime_vf.^theta + (1-delta)*kprime_vf-kprime_vf(ind);

% marginal productivity
margprod = theta*kprime_vf.^(theta-1)+1-delta;

% implied c1
c1_imp = (beta*margprod).^(-1/sigma) .* c2;

%error
EEerror_disc=(c1 - c1_imp)./c1_imp * 100; %percentage deviation
maxEEerror_disc = max(abs(EEerror_disc));


% calculate the relative risk aversion to see the relationship between curvature and EE error.
h = kgrid(2) - kgrid(1); % Calculate the spacing between points
function_derivative = zeros(N,1);
function_derivative(2:N-1) = (kprime_vf(3:N) - kprime_vf(1:N-2)) ./ (2 * h); % Calculate the first derivative
function_derivative(1) = (kprime_vf(2)-kprime_vf(1))/h;
function_derivative(N) = (kprime_vf(N)-kprime_vf(N-1))/h;

second_derivative(2:N-1) = (function_derivative(3:N) - function_derivative(1:N-2)) ./ (2*h); % Calculate the second derivative
second_derivative(1) = (function_derivative(2)-function_derivative(1))/h;
second_derivative(N) = (function_derivative(N)-function_derivative(N-1))/h;

curvature = -(kgrid.*second_derivative)' ./function_derivative;




%Plot for error
    figure_num = figure_num+1;
    figure(figure_num)
    hold on
        yyaxis left
        plot(kgrid,EEerror_disc,'LineWidth',lwidth)
        grid on
        xlabel('Stock of capital')
        ylabel('% error in Euler Equation')
        set(gca,'FontSize',fsize)
    hold on
        yyaxis right
        plot(kgrid, curvature)
        ylabel('Relative Risk Aversion')

%%
% ==============
% 3. Value function iteration with interpolation
% ==============
% set options for fminsearch
options=optimset('MaxIter',5000,'MaxFunEval',5000,'TolFun',1e-12);
% initial guesses
dV=1;
V=zeros(N,1); %zeros(N,1);
iter=0;
% kprime_vf_cont=kprime_vf; %initial guess

        alpha1 = (3-sqrt(5))/2;
        alpha2 = (sqrt(5)-1)/2;
tic

% fminsearch
% while dV>criter_V
%     iter=iter+1;
%     kprimelow=min(kgrid); 
%     kprimehigh=1.3*min(kgrid);
%     for i=1:N % loop over capital today
%         % find maximum over capital tomorrow - now using interpolation
%         [kprime_vf_cont(i),Vnew(i,1)]=fminsearch(@(x) Valuefun(x,kgrid,kgrid(i),theta,sigma,V,delta,beta),kprime_vf_cont(i),options); %%PROBLEM??
%     end
%     Vnew=-Vnew; % take negative as Valuefun is for minimisation
%     % calculate convergence criterion
%     dV=max(max(abs(Vnew-V)));
%     % updated value function,
%     V=Vnew;
%     %disp('dV')
%     %disp(dV)%
%     if iter == 900
%         disp('900');
%     end
% end
t_fmins=toc;


% with Golden Search
dV=1;
ctemp=kgrid.^theta+(1-delta)*kgrid;
V=zeros(N,1);%(ctemp.^(1-sigma)-1)/(1-sigma)
iter=0;
tic

Valuefunpos = @(x)(0);

while dV>criter_V
    iter=iter+1;

    for i=1:N % loop over capital today

        kprimelow=(1-delta) * kgrid(i); 
        kprimehigh=(1-delta)*kgrid(i)+kgrid(i)^theta;
        
        b1= kprimelow+alpha1*(kprimehigh-kprimelow);
        b2= kprimelow+alpha2*(kprimehigh-kprimelow);
        Vlow= ( (kgrid(i)^theta + (1-delta)* kgrid(i)-kprimelow)^( 1-sigma )-1 )/(1-sigma) +...
            beta *Valuefunpos(kprimelow);
        Vhigh=( (kgrid(i)^theta + (1-delta)* kgrid(i)-kprimehigh)^( 1-sigma )-1 )/(1-sigma) +...
            beta *Valuefunpos(kprimehigh);
        Vb1= ( (kgrid(i)^theta + (1-delta)* kgrid(i)-b1)^( 1-sigma )-1 )/(1-sigma) +...
            beta *Valuefunpos(b1);
        Vb2= ( (kgrid(i)^theta + (1-delta)* kgrid(i)-b2)^( 1-sigma )-1 )/(1-sigma) +...
            beta *Valuefunpos(b2);
        dk=1;
        criter_k=1e-12;
        while dk>criter_k
        % use golden search
            if Vb2>Vb1                                            
            kprimelow=b1;  Vlow = Vb1;
            b1=b2;
            Vb1=Vb2;    
            b2=kprimelow+alpha2*(kprimehigh-kprimelow);
            Vb2= ( (kgrid(i)^theta + (1-delta)* kgrid(i)-b2)^( 1-sigma )-1 )/(1-sigma) +...
            beta *Valuefunpos(b2);
            else
            kprimehigh=b2;Vhigh=Vb2;
            b2=b1;Vb2=Vb1;
            b1= kprimelow+alpha1*(kprimehigh-kprimelow);
            Vb1= ( (kgrid(i)^theta + (1-delta)* kgrid(i)-b1)^( 1-sigma )-1 )/(1-sigma) +...
            beta *Valuefunpos(b1);
            end
            dk=abs(Vhigh-Vlow);
        end
        kprime_vf_contG(i)=1/2*(kprimelow+kprimehigh);
        VnewG(i)=1/2*(Vlow+Vhigh);
        
        
        % disp([Vnew(i,1),kprime_VFI_cont(i)])
    end
    
    
    % calculate convergence criterion
    dV=max(max(abs(VnewG-V)));
    % updated value function
    V=VnewG;    
%     disp('dV')
%     disp(dV)
    
    % interpolate
    Valuefunpos = griddedInterpolant(kgrid,V,'spline');
end
toc;
display(iter)

% Plot policy function
figure_num = figure_num +1;
figure(figure_num)
    plot(kgrid,kprime_vf_contG)
    xlabel('Current stock of capital')
    ylabel('Next period''s stock of capital')
    set(gca,'FontSize',fsize)
    hold on
    plot(kgrid, kprime_vf);
    hold off
    
    
    
% Euler equation errors in percent
% consumption vector today
% c1=kgrid.^theta + (1-delta)*kgrid - kprime_vf_contG;
% consumption vector at choice kprime tomorrow
% In two lines, for clarity
        % (1) Get k''(k') by interpolation
        c1=(kgrid.^theta + (1-delta)*kgrid) - kprime_vf_contG;
        kprime2 = interp1(kgrid,kprime_vf_contG ,kprime_vf_contG, 'linear');
        % (2) Get c2
        c2 = kprime_vf_contG.^theta + (1-delta)*kprime_vf_contG - kprime2;

% marginal productivity
margprod=theta*kprime_vf_contG.^(theta-1)+1-delta;

% Implied c1
c1_imp = (beta*margprod).^(-1/sigma) .* c2;

% Error
EEerror_cont = (c1 - c1_imp)./c1_imp * 100;
maxEEerror_cont=max(abs(EEerror_cont));


% calculate the relative risk aversion to see the relationship between curvature and EE error.
h = kgrid(2) - kgrid(1); % Calculate the spacing between points
function_derivative = zeros(N,1);
function_derivative(2:N-1) = (kprime_vf_contG(3:N) - kprime_vf_contG(1:N-2)) ./ (2 * h); % Calculate the first derivative
function_derivative(1) = (kprime_vf_contG(2)-kprime_vf_contG(1))/h;
function_derivative(N) = (kprime_vf_contG(N)-kprime_vf_contG(N-1))/h;

second_derivative(2:N-1) = (function_derivative(3:N) - function_derivative(1:N-2)) ./ (2*h); % Calculate the second derivative
second_derivative(1) = (function_derivative(2)-function_derivative(1))/h;
second_derivative(N) = (function_derivative(N)-function_derivative(N-1))/h;

curvature_linintl = -(kgrid.*second_derivative)' ./function_derivative;

figure_num = figure_num+1;
    figure(figure_num)
    hold on
        yyaxis left
        plot(kgrid,EEerror_cont,'LineWidth',lwidth)
        grid on
        xlabel('Stock of capital')
        ylabel('% error in Euler Equation')
        set(gca,'FontSize',fsize)
    hold on
        yyaxis right
        plot(kgrid, curvature_linintl)
        ylabel('Relative Risk Aversion')


%%
%plotting
