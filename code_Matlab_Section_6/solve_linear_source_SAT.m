%% solve_linear_source_SAT
%
% Description: 
%  Function to numerically solve the linear advection equation with a source term. 
%  Inflow boundary conditions. 
%  The FSBP-SAT method is used on a multi-block structure. 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method. 
%
% Author: Jan Glaubitz 
% Date: Dec 02, 2021 
% 
% INPUT: 
%  x_L, x_R :       left and right boundary of the domain 
%  T :              end time 
%  source :         source term on the RHS of equation 
%  I :              number of blocks 
%  approx_space :  	approximation space (poly, trig, exp, cubic) 
%  K :              dimension of the approximation space 
%  points           data points (equid, Lobatto, Halton, random) 
%  x_eval :         points at which the reference solution is evaluated 
%
% OUTPUT: 
%  x :      grid points 
%  u_num :  numerical solution at grid points 
%  u_ref :  reference solution    
%

function [ x, u_num, u_ref ] = solve_linear_SAT( x_L, x_R, T, source, u_init, I, approx_space, K, points, x_eval )

    %% Set-up the method 

    % Data points and the FSBP operator on the reference block [0,1]
    [ x_ref, w_ref ] = compute_QF( 0, 1, approx_space, points, K ); % grid points and weights on the reference block
    N = length(x_ref); % number of data points
    block_width = (x_R-x_L)/I; % block width 
    [ basis_F, dx_basis_F, span_G, m_G ] = generate_span( 0, 1, approx_space, points, K ); % bases of different spaces 
    [D, P, Q] = compute_FSBP( basis_F, dx_basis_F, x_ref, w_ref ); % FSBP operator 
    D = (1/block_width)*D; P = block_width*P; 
    P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 
    
    % Time step size
    dx_min = min(x_ref(2:end)-x_ref(1:end-1)); % minimum distance between any two neighboring grid points 
    dt = 0.1*(dx_min*block_width); % time-step size
    
    % Global grid points 
    x = zeros(N,I); 
    for i=1:I 
        x(:,i) = x_L + (x_R-x_L)*(i-1)/I + x_ref*block_width;
    end

    % initial data 
    u = u_init(x); % solution values on the global grid 
    if x_eval==0 
        x_eval = x; 
    end
    
    % source term and steady state solution 
    if strcmp( source, '2u')
        source_term = @(x,u) 2*u; % source term
        u_steady = @(x) exp(2*x); % steady state solution 
    elseif strcmp( source, '2xu')
        source_term = @(x,u) 2*x.*u; % source term
        u_steady = @(x) exp(x.^2); % steady state solution 
    else 
       error('Desired source term not yet implemented!') 
    end
    
    u_ref = u_steady(x_eval); % reference solution 

    % Weak enforcement of boundary conditions and coupling 
    SAT = zeros(N,1); % initialize SAT 
    sigma = 1; 

    %% Iterate over time with a 3th-order Runge-Kutta until time T is reached 
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

            % Boundary conditions
            if i == 1
                g_L = 1; % inflow boundary condition 
            else 
                g_L = u(N,i-1); % inter-block coupling 
            end 

            % SAT 
            SAT(1) = -sigma*(u(1,i)-g_L); % SAT term to weakly enforce the boundary condition
            
            %% FSBP
            % 1st update step 
            k1 = u(:,i) + dt*( -D*u(:,i) + source_term(x(:,i),u(:,i)) + P_inv*SAT );  
            % 2nd update step 
            k1 = (3/4)*u(:,i) + (1/4)*k1 + (1/4)*dt*( -D*k1 + source_term(x(:,i),k1) + P_inv*SAT );  
            % 3th update step 
            u_num(:,i) = (1/3)*u(:,i) + (2/3)*k1 + (2/3)*dt*( -D*k1 + source_term(x(:,i),k1) + P_inv*SAT );

        end

        u = u_num; % update solution 
        
    end
    
end