% Create a figure
figure('Visible', 'off');
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 18]);
% Create the first subplot (upper plot)
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(kgrid,kprime_vf_contG,'Linewidth',0.9);
title('Deterministic VFI (Linear Interpolation), Policy Function $k_{t+1}(k_{t})$,$N_k=50, \delta =0,1, \sigma = 2$','interpreter','latex','FontSize',12);
xlabel('Capital $k_t$','interpreter','latex','FontSize',16);
ylabel('Capital $k_{t+1}$','interpreter','latex','FontSize',16);
xline(kbar,'--r');
yline(kbar,'--r');
% Plot the 45-degree line, help to check the steady state point
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k^*'; % The label for the special point
text(6, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(kgrid,VnewG,'Linewidth',0.9);
title('Deterministic VFI(Linear Interpolation) ,Value Function, Nonlinear,$N_k=50, \delta = 0.1. \sigma = 2$','interpreter','latex','FontSize',12);
xlabel('Capital $k_t$','interpreter','latex','FontSize',16);
ylabel('Maximized Utility $V_t$','interpreter','latex','FontSize',16);
grid on;
saveas(gcf,'./figures/vfi_linear_interpolation.png');
hold off; 


% Comparison between VFI in discrete space and VFI with linear interpolation
figure('Visible', 'off');
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 15]);

h1 = plot(kgrid, kprime_vf_contG, 'DisplayName', 'VFI Linear Interpolation', 'Linewidth', 0.9);
hold on;

h2 = plot(kgrid, kprime_vf, 'DisplayName', 'VFI Discrete Maximization', 'Linewidth', 0.9);

title('Deterministic VFI, Policy Function $k_{t+1}(k_{t})$,$N_k=50, \delta =0.1, \sigma = 2$', 'interpreter', 'latex', 'FontSize', 12);
xlabel('Capital $k_t$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('Capital $k_{t+1}$', 'interpreter', 'latex', 'FontSize', 16);

xline(kbar, '--r');
yline(kbar, '--r');
% exclued xline and ylien
legend([h1, h2], 'Location', 'best', 'AutoUpdate', 'off');
saveas(gcf,'./figures/vfi_comparison.png');
hold off;

% Euler Equation Error
figure('Visible', 'off');
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 15]);

plot(kgrid, EEerror_cont, 'Linewidth', 0.9);
hold on;
title('Euler Equation Error, $N_k=50, \delta =0.1, \sigma = 2$', 'interpreter', 'latex', 'FontSize', 12);
xlabel('Capital $k_t$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('Euler Equation Error(\%)', 'interpreter', 'latex', 'FontSize', 16);
saveas(gcf,'./figures/vfi_linearinterp_EulerError.png')


% Euler Equation Error and Curvature

figure('Visible', 'off');
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 15]);

    hold on
        yyaxis left
        plot(kgrid,EEerror_disc,'LineWidth',lwidth)
        grid on
        xlabel('Capital $k_t$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('% error in Euler Equation')
        set(gca,'FontSize',12)
    hold on
        yyaxis right
        plot(kgrid, curvature)
        ylabel('Numerical Policy Curvature')
    title('Euler Equation Error and Policy Curvature',  'FontSize', 12)
saveas(gcf,'./figures/vfi_EulerError_Curvature.png')

% Euler Equation Error and Curvature, VFI linear interpolation
figure('Visible', 'off');
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 15]);
hold on
        yyaxis left
        plot(kgrid,EEerror_cont,'LineWidth',lwidth)
        grid on
        xlabel('Capital $k_t$', 'interpreter', 'latex', 'FontSize', 16)
        ylabel('% error in Euler Equation')
        set(gca,'FontSize',12)
    hold on
        yyaxis right
        plot(kgrid, curvature_linintl)
        ylabel('Numerical Policy Curvature')
title('Euler Equation Error and Policy Curvature, Linear Interpolation',  'FontSize', 12)
saveas(gcf,'./figures/vfi_EulerError_Curvature_VFIlin.png')


% Euler Equation Error for log utility and full depreciation

figure('Visible', 'off');
set(gcf, 'unit', 'centimeters', 'position', [10 5 21 15])
    hold on
        plot(kgrid_log,EEerror_log)
        grid on
        xlabel('Capital $k_t$', 'interpreter', 'latex', 'FontSize', 16)
        ylabel('% error in Euler Equation')
        set(gca,'FontSize',12)
title('Euler Equation Error for the case of Log Utility and Full Depreciation',  'FontSize', 12)
saveas(gcf,'./figures/vfi_EulerError_log.png')


