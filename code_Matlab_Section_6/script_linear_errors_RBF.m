%% script_linear_errors
% 
% Description: 
%  Script to numerically solve the linear advection equation and compare errors 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz and Philipp Öffner 
% Last change Date: Aug 16, 2023 
% Not used inside the paper!

%% Setting up the script 
clc, clear 
 
% Parameters of the problem 
x_L = 0; x_R = 1; % domain boundaries 
T = 0.5; % end time  
%u_init = @(x) ( x > 0 & x < 0.5 ).*( exp(8)*exp( -8./( 1 - (4*x-1).^2 ) ) ); % initial data 
u_init = @(x) (sin( 12*(x-0.1)));
% Shared parameters for the SBP-SAT method  
K = 5; % dimension of approximation space 
x_eval = 0; % evaluation points for reference solution

% Prepare error and loop 
II = []; 
error_L2_poly = []; error_L2_exp = []; 
error_max_poly = []; error_max_exp = []; 

II = [20,40,80,160];
for i=1:length(II) 
    
    I = ceil(II(i)); 
    
    %% Solve the problem using a polynomial function space on Lobatto points 
    x_eval = 0; % evaluation points for reference solution
    approx_space = 'cubic'; % approximation space (poly, trig, exp, cubic)  
    points = 'Halton'; % data points (equid, Lobatto, Halton, random) 
    % Solve problem 
    [ x_poly, u_poly, mass, energy, u_ref ] = solve_linear_SAT_inflow_2( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );
    % Compute errors 
    [ x_ref, w_ref ] = compute_QF( 0, 1, approx_space, points, K ); % grid points and weights on the reference block
    x = x_poly; u = u_poly; 
    error_L2_aux = 0; error_max_aux = 0; 
    for i=1:I 
        error_L2_aux = error_L2_aux + dot(w_ref,(u(:,i)-u_ref(:,i)).^2); 
        error_max_aux = max( error_max_aux, norm( u(:,i)-u_ref(:,i) ) );
    end 
    error_L2_aux = sqrt( error_L2_aux*(x_R-x_L)/I );
    error_L2_poly = [error_L2_poly; error_L2_aux ]; 
    error_max_poly = [error_max_poly; error_max_aux ]; 
    
    %% Solve the problem using a trigonometric function space on equidistant points 
    x_eval = 0; % evaluation points for reference solution 
    approx_space = 'cubic'; % approximation space (poly, trig, exp, cubic)  
    points = 'equid'; % data points (equid, Lobatto, Halton, random) 
    % Solve problem 
    [ x_exp, u_exp, mass, energy, u_ref ] = solve_linear_SAT_inflow_2( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval ); 
    % Compute errors 
    [ x_ref, w_ref ] = compute_QF( 0, 1, approx_space, points, K ); % grid points and weights on the reference block
    x = x_exp; u = u_exp; 
    error_L2_aux = 0; error_max_aux = 0; 
    for i=1:I 
        error_L2_aux = error_L2_aux + dot(w_ref,(u(:,i)-u_ref(:,i)).^2); 
        error_max_aux = max( error_max_aux, norm( u(:,i)-u_ref(:,i) ) );
    end 
    error_L2_aux = sqrt( error_L2_aux*(x_R-x_L)/I );
    error_L2_exp = [error_L2_exp; error_L2_aux ]; 
    error_max_exp = [error_max_exp; error_max_aux ]; 
   
end
    
%% Plot the solutions 

% L2 erros vs I 
figure(1) 
p = plot( II,error_L2_poly,'b^--', II,error_L2_exp,'ro-' ); 
set(p, 'LineWidth',2, 'markersize',12)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
%xlim([40, 81 ]) 
%ylim([ 10^(-15), 1])
xlabel('$I$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_2$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') 
% diagonal line to display rate of convergence 
rate = 2; x = [II(1),II(end)]; y = 1000*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = 'k'; rate_line.LineStyle = ':'; rate_line.LineWidth = 3; 
% diagonal line to display rate of convergence 
rate = 3; x = [II(1),II(end)]; y = 1000*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = [0 0.75 0]; rate_line.LineStyle = '-.'; rate_line.LineWidth = 3; 
lgnd = legend('poly','exp','$2$nd order','$3$th order','Location','best'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 

% Max erros vs I 
figure(2) 
p = plot( II,error_max_poly,'b^--', II,error_max_exp,'ro-' ); 
set(p, 'LineWidth',2, 'markersize',12)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([40, 81 ]) 
%ylim([ 10^(-3), 1])
xlabel('$I$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_\infty$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% diagonal line to display rate of convergence 
rate = 2; x = [II(1),II(end)]; y = 1000*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = 'k'; rate_line.LineStyle = ':'; rate_line.LineWidth = 3; 
% diagonal line to display rate of convergence 
rate = 3; x = [II(1),II(end)]; y = 1000*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = [0 0.75 0]; rate_line.LineStyle = '-.'; rate_line.LineWidth = 3; 
lgnd = legend('poly','exp','$2$nd order','$3$th order','Location','best'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 
