%% script_linear
%
% Description: 
%  Script to numerically solve the linear advection equation 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz and Philipp Öffner 
% Last change Date: Aug 17, 2023
% Used for Test Advection-Diffusion

%% Setting up the script 
clc, clear 
 
%% Parameters of the problem 
x_L = 0; x_R = 0.5; % domain boundaries 
T = 2.0; % end time  
u_init = @(x) (2*x );
u_steady_deriv = @(x,kappa) (1/kappa*(exp(x/kappa))/(exp(1/(2*kappa))-1));
u_steady = @(x,kappa) ((exp(x/kappa)-1)/(exp(1/(2*kappa))-1));
%% Shared parameters for the SBP-SAT method 
K = 5; % dimension of approximation space 
I = 1; % number of blocks  
x_eval = 0; % evaluation points for reference solution
x_eval = linspace(x_L,x_R,1000);

aa=1.0;
kappa=0.2;
% Solve the problem using a cubic function space on equidistant points

approx_space = 'cubic'; % approximation space (cubic, gauss, multi)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_cubic, u_cubic, energy_cubic, u_ref ] = solve_linear_diffusion_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval, kappa, aa);

K = 5;
approx_space = 'gauss'; % approximation space (cubic, gauss, multi)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
%% Solve problem 
[ x_gauss, u_gauss, energy_gauss, u_ref ] = solve_linear_diffusion_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval, kappa, aa );




%% Plots 

% Plot the solutions 
figure(1) 
p = plot( x_eval, u_ref,'k:', x_cubic(:), u_cubic(:),'b--', x_gauss(:), u_gauss(:),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
ylim([-0.01,1.]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','Cubic','Gauss');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% % Plot the mass 
% t = mass_cubic(:,1); % list of times 
% syms x; aux = int(u_init(x),x_L,x_R)
% %aux = 0;
% mass_ref = aux*t.^0; % mass of exact solution over time 
% figure(2) 
% p = plot( t, mass_ref,'k:', t, mass_cubic(:,2),'b-+');%, mass_gauss(:,1),  mass_gauss(:,2),'r-.' ); 
% set(p, 'LineWidth',3)
% set(gca, 'FontSize', 24)  % Increasing ticks fontsize
% %xlim([x(1),x(end)]) 
% %ylim([-1.75,1.75]) 
% xlabel('$t$','Interpreter','latex') 
% ylabel('$\int u \mathrm{d}x$','Interpreter','latex')
% grid on 
% lgnd = legend(p, 'ref','Cubic','Gauss');
% set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')
% 
% Plot the energy 
% t = energy_cubic(:,1); % list of times 
% syms x; aux =int(u_steady(x,kappa)^2,x_L,x_R) +int(2*kappa*(u_steady_deriv(x,kappa)^2),x_L,x_R);
% %aux = 5/8;
% energy_ref = aux*t.^0; % energy of exact solution over time 
% figure(3) 
% p = plot( t, energy_ref,'k:', t, energy_cubic(:,2),'b-+', energy_gauss(:,1), energy_gauss(:,2),'r-.' ); 
% set(p, 'LineWidth',3)
% set(gca, 'FontSize', 24)  % Increasing ticks fontsize
% %xlim([x(1),x(end)]) 
% %ylim([-1.75,1.75]) 
% xlabel('$t$','Interpreter','latex') 
% ylabel('$\int u^2 \mathrm{d}x$','Interpreter','latex')
% grid on 
% lgnd = legend(p, 'ref','Cubic','Gauss');
% set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')