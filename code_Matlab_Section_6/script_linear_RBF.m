%% script_linear
%
% Description: 
%  Script to numerically solve the linear advection equation 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz and Philipp Öffner 
% Last change Date: Aug 16, 2023 
% Used inside the paper in section 6.1!

%% Setting up the script 
clc, clear 
 
%% Parameters of the problem 
x_L = -1; x_R = 1; % domain boundaries 
T = 2.0; % end time  
%u_init = @(x) cos(4*pi*x) + 0.5*sin(40*pi*x); % initial data 
u_init = @(x) exp(-20*x.^2);

%% Shared parameters for the SBP-SAT method 
K = 5; % dimension of approximation space 
I = 20; % number of blocks  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution



% Solve the problem using a cubic function space on equidistant points

approx_space = 'cubic'; % approximation space (cubic, gauss, multi)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_cubic, u_cubic, mass_cubic, energy_cubic, u_ref ] = solve_linear_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );

% 
  approx_space = 'gauss'; % approximation space (cubic, gauss, multi)  
  points = 'equid'; % data points (equid, Lobatto, Halton, random) 
  %% Solve problem 
  [ x_gauss, u_gauss, mass_gauss, energy_gauss, u_ref ] = solve_linear_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );
% % 
% 


%% Plots 

% Plot the solutions 
figure(1) 
p = plot( x_eval, u_ref,'k:',  x_gauss(:), u_gauss(:),'b--', x_cubic(:), u_cubic(:),'r-.' ); 
set(p, 'LineWidth',2)
set(gca, 'FontSize', 12)  % Increasing ticks fontsize
%lim([x(1),x(end)]) 
ylim([-0.1,1.2]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend('ref',    'Gauss','Cubic' );
set(lgnd, 'Interpreter','latex', 'FontSize',12, 'color','none', 'Location','best')

% % Plot the mass % not used inside the manuscipt
t = mass_cubic(:,1); % list of times 
syms x; aux = int(u_init(x),x_L,x_R)
%aux = 0;
mass_ref = aux*t.^0; % mass of exact solution over time 
figure(2) 
p = plot( t, mass_ref,'k:', mass_gauss(:,1),  mass_gauss(:,2),'r--', t, mass_cubic(:,2), 'b-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
ylim([0,0.5]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','Gauss', 'Cubic');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% Plot the energy 
t = energy_cubic(:,1); % list of times 
syms x; aux = int(u_init(x)^2,x_L,x_R)
%aux = 5/8;
energy_ref = aux*t.^0; % energy of exact solution over time 
figure(3) 
p = plot( t, energy_ref,'k:',energy_gauss(:,1), energy_gauss(:,2),'b-.', energy_cubic(:,1), energy_cubic(:,2),'r-+'); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 12)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
ylim([0.2,0.6]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$||u||_2^2$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','Gauss', 'Cubic');
set(lgnd, 'Interpreter','latex', 'FontSize',12, 'color','none', 'Location','best')