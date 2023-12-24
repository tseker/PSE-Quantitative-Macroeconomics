%%
% ==============
% Figures
% ==============

%% Simulation
figure('Name','Simulated time series', 'Visible', 'off')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 18]);
subplot(3,1,1)
hold on
title('Simulation of Capital, Log Linearization')
for sim=1:1
    plot(k_sim_lin(sim,1:T),'Linewidth',1)
end
yline(kbar,'--','Color','red');
text(0.5, kbar+0.01, "$\bar{k}$",'Interpreter','latex')


subplot(3,1,2)
hold on
title('Simulation of Consumption, Log Linearization')
for sim=1:1
    plot(c_sim_lin(sim,1:T),'Linewidth',1)
end
yline(cbar,'--','Color','red');
text(0.5, cbar+0.001, "$\bar{c}$",'Interpreter','latex')



subplot(3,1,3)
hold on
title('Simulation of Labor, Log Linearization')
for sim=1:1
plot(L_sim_lin(sim,1:T),'Linewidth',1)
end
yline(Lbar,'--','Color','red');
text(0.5, Lbar+0.001, "$\bar{L}$",'Interpreter','latex')

saveas(gcf,".\figures\RBC_simulatoin_log_linear.png")

%% Simulation, discrete Markov Chain
figure('Name','Simulated time series', 'Visible', 'off')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 18]);
subplot(3,1,1)
hold on
title('Simulation of Capital, Discrete Markov Process')
for sim=1:1
    plot(k_sim_lin_discrete(sim,1:T),'Linewidth',1)
end
yline(kbar,'--','Color','red');
text(0.5, kbar+0.01, "$\bar{k}$",'Interpreter','latex')


subplot(3,1,2)
hold on
title('Simulation of Consumption, Discrete Markov Process')
for sim=1:1
    plot(c_sim_lin_discrete(sim,1:T),'Linewidth',1)
end
yline(cbar,'--','Color','red');
text(0.5, cbar+0.001, "$\bar{c}$",'Interpreter','latex')



subplot(3,1,3)
hold on
title('Simulation of Labor, Discrete Markov Process')
for sim=1:1
plot(L_sim_lin_discrete(sim,1:T),'Linewidth',1)
end
yline(Lbar,'--','Color','red');
text(0.5, Lbar+0.001, "$\bar{L}$",'Interpreter','latex')

saveas(gcf,".\figures\RBC_simulatoin_log_linear_discrete_Markov.png")





%% Impulse response
figure('Name','Impulse Response', 'Visible', 'off')
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 18]);

subplot(3,3,1)
plot(k_trans);
title('Capital');
xlabel('Period');
yline(kbar,'--','Color','red');
text(0.5, kbar+0.001, "$\bar{K}$",'Interpreter','latex');
hold on;

subplot(3,3,2);
plot(c_trans);
title('Consumption');
xlabel('Period');
yline(cbar,'--','Color','red');
text(0.5, cbar+0.0001, "$\bar{C}$",'Interpreter','latex');
hold on;

subplot(3,3,3)
plot(n_trans);
title('Labor Supply');
xlabel('Period');
yline(Lbar,'--','Color','red');
text(0.5, Lbar+0.00001, "$\bar{N}$",'Interpreter','latex');
hold on;

subplot(3,3,4)
plot(y_trans);
title('Output');
xlabel('Period');
yline(ybar,'--','Color','red');
text(0.5,ybar, "$\bar{Y}$",'Interpreter','latex');
hold on;


subplot(3,3,5)
plot(i_trans);
title('Investment');
xlabel('Period');
yline(ybar-cbar,'--','Color','red');
text(0.5, ybar-cbar+0.00001, "$\bar{I}$",'Interpreter','latex');
hold on;


subplot(3,3,6)
plot(r_trans);
title('Interest Rate');
xlabel('Period');
yline(Rbar,'--','Color','red');
text(0.5, Rbar, "$\bar{r}$",'Interpreter','latex');
hold on;

subplot(3,3,7)
plot(w_trans);
title('Wage');
xlabel('Period');
yline(wbar,'--','Color','red');
text(0.5, wbar, "$\bar{w}$",'Interpreter','latex');
hold on;

subplot(3,3,8)
plot(z_trans);
title('Technology');
xlabel('Period');
yline(1,'--','Color','red');
text(0.5, 1, "$\bar{A}$",'Interpreter','latex');
hold on;

saveas(gcf,".\figures\impulse_response.png");






%% Moments
disp(['Discrete Z Shock Sd.           ','Continuous Z Shock Sd.        ','Real Z Schock Sd.         '])
disp([mean(std(log(Z_sim)')),mean(std(z_cont')),(sigmaepsilon^2/(1-rho^2))^0.5])

disp(['Output Sd.     ','Labor Sd.     ','Consumption Sd.   ','Investment Sd.   '])
% generate the standard deviations
disp([mean(std( y_sim_lin(:,100:T)')),...
    mean(std( L_sim_lin(:,100:T)')),...
    mean(std( c_sim_lin(:,100:T)')),...
    mean(std( inv_sim_lin(:,100:T)')) ])

disp(['Output CV.     ','Labor CV.     ','Consumption CV.   ','Investment CV.   '])
disp([mean(std( y_sim_lin(:,100:T)')./mean(y_sim_lin(:,100:T)')),...
    mean(std( L_sim_lin(:,100:T)')./mean(L_sim_lin(:,100:T)')),...
    mean(std( c_sim_lin(:,100:T)')./mean(c_sim_lin(:,100:T)')),...
    mean(std( inv_sim_lin(:,100:T)')./ mean(inv_sim_lin(:,100:T)') )])


% correlation
rho_emp = zeros(N_sim,1);
rho_emp_cont = zeros(N_sim,1);
corr_Cy_lin = zeros(N_sim,1);
corr_invy_lin = zeros(N_sim,1);
corr_Ly_lin = zeros(N_sim,1);
corr_IC_lin = zeros(N_sim,1);
for i=1:N_sim
    rho_emp(i)=corr(log(Z_sim(i,100:end-1))',log(Z_sim(i,101:end))');
    rho_emp_cont(i)=corr(z_cont(i,100:end-1)',z_cont(i,101:end)');
    corr_Cy_lin(i)=corr(y_sim_lin(i,100:T)',c_sim_lin(i,100:T)');
    corr_invy_lin(i)=corr(y_sim_lin(i,100:T)',inv_sim_lin(i,100:T)');
    corr_Ly_lin(i)=corr(y_sim_lin(i,100:T)',L_sim_lin(i,100:T)');
    corr_IC_lin(i)=corr(inv_sim_lin(i,100:T)',c_sim_lin(i,100:T)');
end


% generate the mean of the moments
rho_emp_mean = mean(rho_emp);
rho_emp_cont_mean = mean(rho_emp_cont);
corr_Cy_lin_mean = mean(corr_Cy_lin);
corr_invy_lin_mean = mean(corr_invy_lin);
corr_Ly_lin_mean = mean(corr_Ly_lin);
corr_IC_lin_mean = mean(corr_IC_lin);

disp(['Output-Consumption Corr.     ','Output-Investment Corr.     ','Output-Labor Corr.    ', 'Investment-Consumption Corr.    '])
% generate the standard deviations
disp([corr_Cy_lin_mean,corr_invy_lin_mean,corr_Ly_lin_mean, corr_IC_lin_mean])

%% Moments, discrete
disp(['Discrete Z Shock Sd.           ','Continuous Z Shock Sd.        ','Real Z Schock Sd.         '])
disp([mean(std(log(Z_sim)')),mean(std(z_cont')),(sigmaepsilon^2/(1-rho^2))^0.5])

disp(['Output Sd.     ','Labor Sd.     ','Consumption Sd.   ','Investment Sd.   '])
% generate the standard deviations
disp([mean(std( y_sim_lin_discrete(:,100:T)')),...
    mean(std(L_sim_lin_discrete(:,100:T)')),...
    mean(std(c_sim_lin_discrete(:,100:T)')),...
    mean(std(inv_sim_lin_discrete(:,100:T)'))])


disp(['Output CV.     ','Labor CV.     ','Consumption CV.   ','Investment CV.   '])
disp([mean(std( y_sim_lin_discrete(:,100:T)')./mean(y_sim_lin_discrete(:,100:T)')),...
    mean(std( L_sim_lin_discrete(:,100:T)')./mean(L_sim_lin_discrete(:,100:T)')),...
    mean(std( c_sim_lin_discrete(:,100:T)')./mean(c_sim_lin_discrete(:,100:T)')),...
    mean(std( inv_sim_lin_discrete(:,100:T)')./ mean(inv_sim_lin_discrete(:,100:T)') )])



% correlation
rho_emp = zeros(N_sim,1);
rho_emp_cont = zeros(N_sim,1);
corr_Cy_lin_disc = zeros(N_sim,1);
corr_invy_lin_disc = zeros(N_sim,1);
corr_Ly_lin_disc = zeros(N_sim,1);
for i=1:N_sim
    rho_emp(i)=corr(log(Z_sim(i,100:end-1))',log(Z_sim(i,101:end))');
    rho_emp_cont(i)=corr(z_cont(i,100:end-1)',z_cont(i,101:end)');
    corr_Cy_lin_disc(i)=corr(y_sim_lin_discrete(i,100:T)',c_sim_lin_discrete(i,100:T)');
    corr_invy_lin_disc(i)=corr(y_sim_lin_discrete(i,100:T)',inv_sim_lin_discrete(i,100:T)');
    corr_Ly_lin_disc(i)=corr(y_sim_lin_discrete(i,100:T)',L_sim_lin_discrete(i,100:T)');
end


% generate the mean of the moments
rho_emp_mean = mean(rho_emp);
rho_emp_cont_mean = mean(rho_emp_cont);
corr_Cy_lin_mean = mean(corr_Cy_lin_disc);
corr_invy_lin_mean = mean(corr_invy_lin_disc);
corr_Ly_lin_mean = mean(corr_Ly_lin_disc);

disp(['Output-Consumption Corr.     ','Output-Investment Corr.     ','Output-Labor Corr.    '])
% generate the standard deviations
disp([corr_Cy_lin_mean,corr_invy_lin_mean,corr_Ly_lin_mean])



