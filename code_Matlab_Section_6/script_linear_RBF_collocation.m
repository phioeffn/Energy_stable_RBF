%% script_linear
%
% Description: 
%  Script to numerically solve the linear advection equation 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz and Philipp Öffner 
% Last change date: Aug 16, 2023 
% Used inside the paper to calculate the first example!

%% Setting up the script 
clc, clear 
 
%%%% For the classical stuff

% Script to compare collocation and FSBP operators The first part is from
% the Paper by Glaubitz and Gelb. The Second one is from our paper here. 

% We use the domain [-1,1] 
clear, clc, close all 

%% Setting up common variables 
Init_C = 'exp'; % sin, exp, disc
BC = 'periodic'; % inflow, periodic
T = 2; % final time 
basis = 'cubic'; % G, MQ, IQ, cubic, quintic
ep = 5; % shape parameter
N = 15; % number of points 
d = -1; % polynomial degree 
points = 'equid'; % equid, random
CFL = 0.1; % CFL number 
integration = 'exact'; % way integration is performed (exact, trapez, Gauss)

x = linspace(-1,1,N)';
xxx = sort(x,'ascend');
%% set up RBF and IC 
rbf = basis_function( basis );
IC = initial_cond( Init_C ); 
u0 = IC(x);

%% routine for strong and weak RBF method 
[u_strong, m_strong, e_strong] = linear_strong_RBF( BC, T, CFL, x, IC, rbf, ep); % strong RBF




%% Parameters of the problem 
x_L = -1; x_R = 1; % domain boundaries 
T = 2.0; % end time  
%u_init = @(x) cos(4*pi*x) + 0.5*sin(40*pi*x); % initial data 
u_init = @(x) exp(-20*x.^2);

%% Shared parameters for the SBP-SAT method 
K = 14; % dimension of approximation space 
I = 1; % number of blocks  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution



% Solve the problem using a cubic function space on equidistant points

approx_space = 'cubic'; % approximation space (cubic, gauss, multi)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_cubic, u_cubic, mass_cubic, energy_cubic, u_ref ] = solve_linear_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );



%% Plots 

% Plot the solutions 
figure(1) 
p = plot( x_eval, u_ref,'k:',xxx(:),u_strong(:), 'r--', x_cubic(:), u_cubic(:),'b-.');%,x_cubic(:), u_cubic(:),'r--'); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 12)  % Increasing ticks fontsize
%lim([x(1),x(end)]) 
ylim([-0.1,1.2]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend('ref', 'Col.','SBP');
set(lgnd, 'Interpreter','latex', 'FontSize',12, 'color','none', 'Location','best')


% Plot the energy 
t = energy_cubic(:,1); % list of times 
syms x; aux = int(u_init(x)^2,x_L,x_R)
%aux = 5/8;
energy_ref = aux*t.^0; % energy of exact solution over time 
figure(3) 
p = plot( t, energy_ref,'k:',e_strong(:,1), e_strong(:,2),'r-.', energy_cubic(:,1), energy_cubic(:,2),'b-+'); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 12)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([0.2802,0.2803]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$||u||_2^2$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','Col.', 'SBP');
set(lgnd, 'Interpreter','latex', 'FontSize',12, 'color','none', 'Location','best')