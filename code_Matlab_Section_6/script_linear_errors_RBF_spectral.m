%% script_linear_errors
%
% Description: 
%  Script to numerically solve the linear advection equation and compare errors 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz and Philipp Öffner
% Date: Aug 16, 2023
% USed for the spectral convergence investigation 

%% Setting up the script 
clc, clear 
 
% Parameters of the problem 
x_L = 0; x_R = 1.; % domain boundaries 
T = 0.5; % end time  
u_init = @(x) ( x > 0 & x < 0.5 ).*( exp(8)*exp( -8./( 1 - (4*x-1).^2 ) ) ); % initial data 


% Shared parameters for the SBP-SAT method  
I= 20; % Number of blocks 
x_eval = 0; % evaluation points for reference solution

% Prepare error and loop 
k = []; 
error_L2_poly = []; error_L2_exp = []; 
error_max_poly = []; error_max_exp = []; 

kk= [4,5,6,7];
for i=1:length(kk) 
    
    K = ceil(kk(i)); 
    
    %% Solve the problem using a cubiv function space on equidistant points 
    x_eval = 0; % evaluation points for reference solution
    approx_space = 'cubic'; % approximation space (poly, trig, exp, cubic)  
    points = 'equid'; % data points (equid, Lobatto, Halton, random) 
    % Solve problem 
    [ x_poly, u_poly, mass, energy, u_ref ] = solve_linear_SAT_inflow( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );
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
    

   
end
    
%% Plot the solutions 

% L2 erros vs I 
figure(1) 
%p = plot( kk,log(error_L2_poly),'b^--');%,kk,log(error_L2_exp),'ro-' );
p=loglog(kk,error_L2_poly,'b^--')
set(p, 'LineWidth',2, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
%xlim([40, 81 ]) 
%ylim([ 10^(-15), 1])
xlabel('$K$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_2$','Interpreter','latex')

 grid on 

% Max erros vs I 
figure(2) 
%p = plot( kk, log(error_max_poly),'b^--');%,

p=loglog(kk,error_max_poly,'b^--');%(, "ro-");%kk,log(error_max_exp),'ro-' ); 
set(p, 'LineWidth',2, 'markersize',10)
set(gca, 'FontSize',20)  % Increasing ticks fontsize
%xlim([40, 81 ]) 
xlabel('$K$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_\infty$','Interpreter','latex')

grid on 
