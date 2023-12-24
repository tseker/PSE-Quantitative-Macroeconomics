
%% Broyden, Infinitely living

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(x(1:T,1),'Linewidth',2);
title('NCGM Path of Capital Stock with Broyden Method','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k_*'; % The label for the special point
line([min(x), 150], [specialPointK, specialPointK], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(x(T+1:2*T,1),'Color', 'black','LineWidth',2);
title('NCGM Path of Consumption with with Broyden Method','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
specialPointC = cbar; % Replace with the y-coordinate of the special point
specialPointLabelC = 'c_*'; % The label for the special point
line([min(x), 150], [specialPointC, specialPointC], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointC, specialPointLabelC, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
saveas(gcf,'broyden.png');
hold off; 


asdfasdf




%% Multiple Shooting, Infinitely Living
% Plotting
figure;
% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_path_ms,'Linewidth',2);
title('NCGM Path of Capital Stock with Multiple Shooting','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k_*'; % The label for the special point
line([min(x), 150], [specialPointK, specialPointK], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_path_ms, 'Color', 'black','LineWidth',2);
title('NCGM Path of Consumption with Multiple Shooting','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
specialPointC = cbar; % Replace with the y-coordinate of the special point
specialPointLabelC = 'c_*'; % The label for the special point
line([min(x), 150], [specialPointC, specialPointC], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointC, specialPointLabelC, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
saveas(gcf,'mshooting.png');
hold off; 







%% finitely Living

% ==============
% fsolve plotting
% ==============

% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_finite_living_T10,'Linewidth',2);
title('fsolve: Deterministic Path of Capital Stock when T = 10','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;



% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_finite_living_T10, 'Color', 'black','LineWidth',2);
title('fsolve: Deterministic Path of Consumption when T = 10','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T10.png');
hold off; 


% T = 100 Case
% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(real(k_finite_living_T100),'Linewidth',2);
title('fsolve: Deterministic Path of Capital Stock when T = 100','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;


% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(real(c_finite_living_T100), 'Color', 'black','LineWidth',2);
title('fsolve: Deterministic Path of Consumption when T = 100','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T100.png');
hold off; 

% T = 200 Case
% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(real(k_finite_living_T200),'Linewidth',2);
title('fsolve: Deterministic Path of Capital Stock when T = 200','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;

% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_finite_living_T200, 'Color', 'black','LineWidth',2);
title('fsolve: Deterministic Path of Consumption when T = 200','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T200.png');
hold off; 

% ==============
% Broyden plotting
% ==============

% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_trans_br_T10,'Linewidth',2);
title('Broyden: Deterministic Path of Capital Stock when T = 10','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;



% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_trans_br_T10, 'Color', 'black','LineWidth',2);
title('Broyden: Deterministic Path of Consumption when T = 10','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T10_br.png');
hold off; 


% T = 100 Case
% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(real(k_trans_br_T100),'Linewidth',2);
title('Broyden: Deterministic Path of Capital Stock when T = 100','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;


% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(real(c_trans_br_T100), 'Color', 'black','LineWidth',2);
title('Broyden: Deterministic Path of Consumption when T = 100','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T100_br.png');
hold off; 

% T = 200 Case
% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(real(k_trans_br_T200),'Linewidth',2);
title('Broyden: Deterministic Path of Capital Stock when T = 200','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;

% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_trans_br_T200, 'Color', 'black','LineWidth',2);
title('Broyden: Deterministic Path of Consumption when T = 200','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T200_br.png');
hold off; 

% ==============
% Multiple Shooting
% ==============

% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_path_ms_t10,'Linewidth',2);
title('Multiple Shooting: Deterministic Path of Capital Stock when T = 10','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;



% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_path_ms_t10, 'Color', 'black','LineWidth',2);
title('Multiple Shooting: Deterministic Path of Consumption when T = 10','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T10_ms.png');
hold off; 


% T = 100 Case
% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(real(k_path_ms_t100),'Linewidth',2);
title('Multiple Shooting: Deterministic Path of Capital Stock when T = 100','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;


% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(real(c_path_ms_t100), 'Color', 'black','LineWidth',2);
title('Multiple Shooting: Deterministic Path of Consumption when T = 100','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T100_ms.png');
hold off; 

% T = 200 Case
% Create a figure
figure;

% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(real(k_path_ms_t200),'Linewidth',2);
title('Multiple Shooting: Deterministic Path of Capital Stock when T = 200','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;

% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_path_ms_t200, 'Color', 'black','LineWidth',2);
title('Multiple Shooting: Deterministic Path of Consumption when T = 200','interpreter','latex','FontSize',12);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'finite_T200_ms.png');
hold off; 

