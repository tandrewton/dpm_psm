%% plot the force
close all; clear;
set(0,'DefaultFigureWindowStyle','docked')
epsilon = 1;
a = 0.1;
b = 0.3;
g = @(x,a,b)(2*epsilon*(1-x).*(x < 1 + a) ... 
    - 2*epsilon*a/(b-a)*(1 + b - x).*(x >= 1 + a));
domain = linspace(0.5, 1+b, 1000);
figure(1)
hold on
yline(0,'k--','handlevisibility', 'off','LineWidth', 1)
plot(domain, g(domain,0,b),'green--','linewidth', 3, 'DisplayName', '$~\epsilon_2^*$ = 0.0')
plot(domain, g(domain,0.01,b),'green-','linewidth', 3, 'DisplayName', '$~\epsilon_2^*$ = 0.01')
plot(domain, g(domain,0.015,b),'magenta--','linewidth', 3, 'DisplayName', '$~\epsilon_1^*$ = 0.015')
plot(domain, g(domain,0.035,b),'magenta-','linewidth', 3, 'DisplayName', '$~\epsilon_1^*$ = 0.035')

set(gca, 'FontSize', 24)
xlabel('$d/\sigma$','fontsize', 30, 'interpreter', 'latex')
ylabel('$f_{int}\sigma/u_0$', 'fontsize', 30, 'interpreter', 'latex')
legend('fontsize', 24, 'interpreter', 'latex')

xticks([0.9,1,1.1,1.2,1.3])
xlim([0.85,1.3])
box on
ax = gca;
ax.TickLength = [0.025 0.025];
ax.LineWidth = 1;

%% plot the energy
figure(2); hold on;
g = @(x,a,b)(epsilon*((1-x).^2 - a*b).*(x < 1 + a) ... 
    - epsilon*a/(b-a)*(1 + b - x).^2.*(x >= 1 + a));
yline(0,'k--','handlevisibility', 'off','LineWidth', 1)
plot(domain, g(domain,0,b),'green--','linewidth', 3,'DisplayName', '$~\epsilon_2^*$ = 0.0')
plot(domain, g(domain,0.01,b),'green-','linewidth', 3, 'DisplayName', '$~\epsilon_2^*$ = 0.01')
plot(domain, g(domain,0.015,b),'magenta--','linewidth', 3, 'DisplayName', '$~\epsilon_1^*$ = 0.015')
plot(domain, g(domain,0.035,b),'magenta-','linewidth', 3, 'DisplayName', '$~\epsilon_1^*$ = 0.035')

set(gca, 'FontSize', 24)
xlabel('$d/\sigma$','fontsize', 30, 'interpreter', 'latex')
ylabel('$u_{int}/u_0$', 'fontsize', 30, 'interpreter', 'latex')
legend('fontsize', 24, 'interpreter', 'latex')

xticks([0.9,1,1.1,1.2,1.3])
xlim([0.85,1.3])
ylim([-0.012, 0.025])
yticks([-0.01, 0, 0.01, 0.02])
box on  
ax = gca;
ax.TickLength = [0.025 0.025];
ax.LineWidth = 1;
% Legend
legend('show');
legend('Location', 'northeast');
hold off;