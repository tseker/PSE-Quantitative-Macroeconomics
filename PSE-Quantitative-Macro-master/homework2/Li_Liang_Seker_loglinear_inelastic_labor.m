% ==========================
% Question 2: Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
close all
clear all

% change location to the folder where this m file is saved
mfile_name          = mfilename('fullpath');  % This line of code retrieves the full path of the current MATLAB file and assigns it to the variable mfile_name. The mfilename function with the argument 'fullpath' returns the full path of the currently executing file, including the file name and extension. The full path is then assigned to the variable mfile_name.
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);


% ============
% parameters
% ============
alpha = 0.3; % capital share
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 1.0001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.025;

% ============
% options and convergence criteria
% ============
T=150; % periods for transition

%mean of capital non-stochastic steady state
kbar = (alpha*beta/(1-beta*(1-delta)))^(1/(1-alpha)); %write feasibility constraint at the steady-state and EE. Then, replace c_bar from feasibiltiy constraint into the EE, and obtain k_bar. 
% initial level of capital in the transition
k_0=kbar*0.9;

%%
% =====================
% Log-linearization
% =====================

% some ratios as function of parameters
ybar = kbar^(alpha);
cbar = ybar-delta*kbar;
ckrat=cbar/kbar; %c/k ratio
R=1/beta;
margprod=R-1+delta;

% a. write system as A E[y_t+1]+B y_t=0

% Euler equation
% -gamma*chat_t=(-gamma)*chat_t+1
% +beta*margprod*khat_t+1+beta*margprod*zbar*zhat_t+1

% becomes
% (-gamma)*chat_t-beta*margprod*khat_t+1  -  beta*margprod*zhat_t+1+gamma chat_t=0

% res constraint
% khat_t+1 - ckrat*chat_t - (margprod+(1-delta)) khat_t -kbar^theta*zbar*zhat_t=0

% order c,k,z
A=[-sigma beta*alpha*(alpha-1)*kbar^(alpha-1); 0 1];

B= [-sigma 0; -ckrat alpha*kbar^(alpha-1)+(1-delta)];

D = A\B; %inv(A)*B can be slower;

% note that these are right-hand eigenvectors
[ev, lambda]=eig(D); %Note: [V,D] = eig(A) produces a diagonal matrix D of eigenvalues and 
                            %a full matrix V whose columns are the corresponding eigenvectors  
                            %so that A*V = V*D
aaa=inv(ev);


BKcond=sum(abs(diag(lambda))>1);
%If all the eigenvalues are less than 1 in absolute value then the system is stable. 
% If all the eigenvalues are greater than 1 in absolute value then the system is unstable. 
% If at least one eigenvalue is less than 1 in absolute value the system is saddle-path stable.

if BKcond~=1
    disp('BK conditions not satisfied')
else
    indic=find(abs(diag(lambda))>1);
    indic1=find(abs(diag(lambda))<=1);
    polfunc_temp=aaa(indic,:);
    % this gives the consumption policy function
    polfunc = -polfunc_temp(1,2)/polfunc_temp(1,1); % policy function for chat(t+1)
    % Recall that we found policy function as -p12/p11*khat(t+1)
end


% calculate the deterministic transition
% 
c_lin = zeros(T,1);
k_lin = zeros(T,1);
% k_lin(1)=0.9*kbar; %k0
% for t=2:T
%     c_lin(t-1) = polfunc*((k_lin(t-1)-kbar)/kbar)*cbar+cbar;
%     k_lin(t)=k_lin(t-1)^alpha+(1-delta)*k_lin(t-1)-c_lin(t-1);
% end
% c_lin(T) = polfunc*((k_lin(T)-kbar)/kbar)*cbar+cbar;
k_lin(1)=-0.1; %k0
for t=2:T
    c_lin(t-1) = polfunc*k_lin(t-1);
    k_lin(t)=  -cbar/kbar* c_lin(t-1) + (alpha * kbar^(alpha-1) +(1-delta)) * k_lin(t-1);
end

c_lin(T) = polfunc*k_lin(T);
c_lin = c_lin.*cbar + cbar;
k_lin = k_lin.*kbar + kbar;


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
criter_V = 0.0001;
c_0n=(c_0l+c_0u)/2;
errorn=cT(c_0n,k_0,T,alpha,delta,sigma,beta)-cbar;
while abs(errorn)>criter_V
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


% find the path

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





%%
% ==============
% Figures
% ==============

figure
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_lin,'Linewidth',2);
title('Simulated transition for capital stock','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k_*'; % The label for the special point
line([min(k_lin), 150], [specialPointK, specialPointK], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_lin,'Color', 'black','LineWidth',2);
title('Simulated transition for consumption','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
specialPointC = cbar; % Replace with the y-coordinate of the special point
specialPointLabelC = 'c_*'; % The label for the special point
line([min(c_lin), 150], [specialPointC, specialPointC], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointC, specialPointLabelC, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
saveas(gcf,'.\figures\q2_loglinpath.png');



% Multiple shooting figure

figure
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_path_ms,'Linewidth',2);
title('Simulated transition for capital stock, multiple shooting','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k_*'; % The label for the special point
line([min(k_path_ms), 150], [specialPointK, specialPointK], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_path_ms,'Color', 'black','LineWidth',2);
title('Simulated transition for consumption, multiple shooting','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
specialPointC = cbar; % Replace with the y-coordinate of the special point
specialPointLabelC = 'c_*'; % The label for the special point
line([min(c_path_ms), 150], [specialPointC, specialPointC], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointC, specialPointLabelC, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
saveas(gcf,'.\figures\q2_multipleshooting.png');




