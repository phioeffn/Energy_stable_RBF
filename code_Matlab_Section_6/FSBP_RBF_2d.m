%% script_linear_2d
%
% Description: 
%  Script to numerically solve the linear advection equation in 2d
%  ? boundary conditions 
%  The FSBP-RBF method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz & Philipp Oeffner
% Date: March 04, 2022 

%% Setting up the script 
clc, clear 
 
% Parameters of the problem 
x_L = 0; x_R = 1; % domain boundaries 
T = 0.5; % end time  

% Shared parameters for the SBP-SAT method 
K = 14; % number of center points  %14 is used in the paper 
I = 1; % number of blocks 
approx_space = 'cubic'; % kernel 
points = 'equid'; % type of points
sigma = 1.0; 


%% 1D matrices 

[ x_ref, w_ref ] = compute_QF( 0, 1, approx_space, points, K ); % grid points and weights on the 
n = length(x_ref); % number of data points
block_width = (x_R-x_L)/I; % block width
[ basis_F, dx_basis_F, span_G, m_G ] = generate_span( 0, 1, approx_space, points, K ); % bases of different spaces 
[D, P, Q] = compute_FSBP( basis_F, dx_basis_F, x_ref, w_ref ); % FSBP operator 
D = (1/block_width)*D; 
P = block_width*P; 
P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 

%% 2D matrices 

Dx = kron(D,eye(n)); 
Dy = kron(eye(n),D); 

P_2d = kron(P,P); % Is doing the same
P_2d_inv = inv(P_2d); 

%% 2D grid points 
 
% tensor product points in 2d 
n =length(x_ref); 
x = x_ref; 
[xx, yy] = meshgrid(x, x); 
X = [ reshape(xx, numel(xx),1)'; 
reshape(yy, numel(yy),1)']'; 


%% Initial condition 

%IC = @(x,y) (x>=0.5).*sin(4*pi*x).*(1-0.5*sin(2*pi*y)); % initial condition 
IC = @(x,y) exp(-20*( (x-0.25).^2 + (y-0.25).^2 ) ); % initial condition 
u_init = IC(X(:,1),X(:,2)); % IC for numerical solution 

uu_init = reshape( u_init, n, n ); % reshape into matrix for plotting 
figure(1) 
s = mesh(xx,yy,uu_init); 
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 
%zlim([-2,2]) 


%% Time integration 
% Time step size
dx_min = min(x_ref(2:end)-x_ref(1:end-1)); % minimum distance between any two neighboring grid points 
dt = 0.01*(dx_min*block_width); % time-step size 
u = u_init; 
 mass = []; energy = []; % mass and energy over time 
t = 0; 
while (t<T)  
    
	% time stepping 
   	if T-t<dt 
      	dt = T-t; 
    else
      	t = t+dt
    end
   
    % loop over block 
    for i = 1:I
    
        % SAT 
        SAT_x = zeros(n,n); 
        SAT_y = zeros(n,n); 
        uu = reshape( u, n, n ); % reshape into matrix for plotting 
        g_L = zeros(n,1); % zero inflow 
        g_U = zeros(1,n); % zero inflow
       % g_L = uu(:,n); % for periodic BC 
       % g_U = uu(n,:); % for periodic BC 
        SAT_x(:,1) = -sigma*0.1*(uu(:,1)-g_L); % SAT term to weakly enforce the boundary condition
        SAT_x = P_2d_inv*SAT_x(:);    
        %SAT = P_inv(1,1)*SAT(:); 
        SAT_y(1,:) = -sigma*0.1*(uu(1,:)-g_U); % SAT term to weakly enforce the boundary condition
        SAT_y = P_2d_inv*SAT_y(:); 
        
        % 1st update step 
        k1 = u + dt*( -0.5*Dx*u -Dy*u + SAT_x + SAT_y );  
        % 2nd update step 
        k2 = (3/4)*u + (1/4)*k1 + (1/4)*dt*( -0.5*Dx*k1 -Dy*k1 + SAT_x + SAT_y );  
        % 3th update step 
        u_num = (1/3)*u + (2/3)*k2 + (2/3)*dt*( -0.5*Dx*k2 -Dy*k1 + SAT_x + SAT_y );

    end
        
    u = u_num; % update solution   
    
    % Mass and energyy 
    
    mass_aux = 0; energy_aux = 0; 
     
      mass_aux = dot( ones(n^2,1), P_2d*u(:,i)); % compute mass 
      energy_aux = dot( u(:,i), P_2d*u(:,i) ); % compute energy
   
        % save mass and energy
        mass = [mass; t, mass_aux]; % save mass
        energy = [energy; t, energy_aux]; % save 
end 

uu = reshape( u, n, n ); % reshape into matrix for plotting 
figure(2) 
s = mesh(xx,yy,uu); 
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 
%zlim([-2,2]) 

%u_ana = @(x,y,t) ( x > T ).*(sin(4*pi*(x-T).*(1-0.5*sin(2*pi*y))));
%u_ref = u_ana(X(:,1),X(:,2),T); % IC for numerical solution 

%uu_ref = reshape( u_ref, n, n ); % reshape into matrix for plotting 
%figure(3) 
%s = mesh(xx,yy,uu_ref); 
%s.EdgeColor = 'interp'; 
%set(s, 'LineWidth',1.5); 
%set(gca, 'FontSize', 20)  % Increasing ticks fontsize 
%xlabel('$x$','Interpreter','latex') 
%ylabel('$y$','Interpreter','latex') 
%zlabel('$u$','Interpreter','latex') 
%zlim([-2,2]) 







% % Plot the mass 
t = mass(:,1); % list of times 
syms x; aux = integral2(IC,x_L,x_R,x_L,x_R)
%aux = 0;
mass_ref = aux*t.^0; % mass of exact solution over time 
figure(3) 
p = plot( t, mass_ref,'k:',  mass(:,1), mass(:,2), 'b-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([0.15,0.16]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref', 'Cubic');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

IC_2 = @(x,y) exp(-20*( (x-0.25).^2 + (y-0.25).^2 ) );

% Plot the energy 
t = energy(:,1); % list of times 
syms x; aux = integral2(IC_2,x_L,x_R,x_L,x_R)
%aux = 5/8;
energy_ref = aux*t.^0; % energy of exact solution over time 
figure(4) 
p = plot( t, energy_ref,'k:', energy(:,1), energy(:,2),'b-+'); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlim([0,0.5]) 
ylim([0.07,0.09]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$$||u||_2^2$$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref', 'Cubic');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')












