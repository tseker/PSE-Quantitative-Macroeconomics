%% 
% =====================
% Endogenous Labor Supply Model
% =====================

close all
clear all

% change location to the folder where this m file is saved
mfile_name          = mfilename('fullpath');  % This line of code retrieves the full path of the current MATLAB file and assigns it to the variable mfile_name. The mfilename function with the argument 'fullpath' returns the full path of the currently executing file, including the file name and extension. The full path is then assigned to the variable mfile_name.
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);


% parameters

beta = 0.99;
alpha = 0.3;
sigma = 1.0001;
delta = 0.025;
mu = 1;


%========== Update steady state value of capital, labor and consumption



% 1 = beta*[(1-delta) + alpha * kbar^(alpha-1) * lbar^(1-alpha) ]
% theta * lbar ^(mu+alpha) = beta * (1-alpha) * cbar ^(-sigma) * kbar^(-alpha)
% kbar = (1-delta) * kbar + kbar^alpha * lbar^(1-alpha) - cbar

% Here we denote m = lbar/kbar
% From 1 = beta*[(1-delta) + alpha * kbar^(alpha-1) * lbar^(1-alpha) ], we get:

m = ( alpha/(1/beta-(1-delta)) )^(1/(1-alpha));
lbar = 1/3;

% Calibration set lbar to 1/3, mu = 1
% From the three equations above, we get:
% theta * lbar^(mu+sigma) = beta*(1-alpha)*(m^alpha-delta*m)^(-sigma) * m^alpha

theta = lbar^(-mu-sigma) *(1-alpha)*(m^alpha-delta*m)^(-sigma) * m^alpha;

% Steady State kbar and cbar

kbar = m * lbar ; % recover kbar from m
cbar = -delta* kbar + kbar^alpha * lbar^(1-alpha);



% ============== log linearize =======================

% A x_{t+1} = B x_t, note that here xt is the deviation from the steady state.
% xt [kt;ct;lt]
A = [1, 0, 0;
    beta*alpha*(alpha-1)*lbar^(1-alpha)*kbar^(alpha-1), -sigma, beta*alpha*( 1-alpha ) * kbar^(alpha-1) * lbar^(1-alpha);
    0, 0, 0];

B = [ (1-delta)+ alpha * kbar^(alpha-1) * lbar^(1-alpha), - cbar/kbar, (1-alpha)*kbar^(alpha-1)* lbar^(1-alpha);
       0, -sigma, 0;
       -alpha, sigma, mu + alpha];
   
[polic_func_endog_l, law_motion]  = solab(A,B,1);

k0 = 0.9 * kbar;


%% Time Series
T = 150;
k1 = 0.9 * kbar;
k_lin = zeros(T,1); % k1, k2,..., kT, linearized
c_lin = zeros(T,1); % linearized c
l_lin = zeros(T,1); 


k_lin(1) =( k1-kbar)/kbar;
k(1) = k1;

for t = 1:T
   temp = polic_func_endog_l * k_lin(t);
   c_lin(t) = temp(1);
   l_lin(t) = temp(2);
   if t ~= T
        k_lin(t+1) = law_motion * k_lin(t);
   end
end

% recover k, c, and l form the k, c and l
k = k_lin .* kbar + kbar;
c = c_lin .* cbar + cbar;
l = l_lin .* lbar + lbar;

set(gcf, 'unit', 'centimeters', 'position', [10 5 21 18]);

subplot(3,1,1)
plot(k,'Linewidth',2);
title('Capital stock, endogenous labor supply', 'interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = '$\bar{k}$'; % The label for the special point
line([1, 150], [specialPointK, specialPointK], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointK, specialPointLabelK,'interpreter','latex','FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');





subplot(3,1,2)
plot(c,'o','Linewidth',2,'color',[0.6350 0.0780 0.1840]);
title('Consumption, endogenous labor supply', 'interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
specialPointC = cbar; % Replace with the y-coordinate of the special point
specialPointLabelC = '$\bar{c}$'; % The label for the special point
line([1, 150], [specialPointC, specialPointC], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointC, specialPointLabelC,'interpreter','latex','FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');



subplot(3,1,3)
plot(l,"+",'Linewidth',2, 'color',[0.9290 0.6940 0.1250]);
title('Labor, endogenous labor supply', 'interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Labor $l_t$','interpreter','latex','FontSize',16);

grid on;
specialPointL = lbar; % Replace with the y-coordinate of the special point
specialPointLabelL = '$\bar{l}$'; % The label for the special point
line([1, 150], [specialPointL, specialPointL], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointL, specialPointLabelL,'interpreter','latex','FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

saveas(gcf,'.\figures\Endongenous_Labor_Supply_Transition.png');














