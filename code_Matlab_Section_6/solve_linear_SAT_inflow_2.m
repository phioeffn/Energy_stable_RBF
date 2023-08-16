%% solve_linear_SAT
%
% Description: 
%  Function to solve the linear advection equation with periodic boundary
%  conditions using a single/multi-block FSBP-SAT method. 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: Dec 02, 2021 
% 
% INPUT: 
%  x_L, x_R :       left and right boundary of the domain 
%  T :              end time 
%  u_init :         initial condition 
%  I :              number of blocks 
%  approx_space :  	approximation space (poly, trig, exp, cubic) 
%  K :              dimension of the approximation space 
%  points           data points (equid, Lobatto, Halton, random) 
%  x_eval :         points at which the reference solution is evaluated 
%
%
% OUTPUT: 
%  x :              grid points 
%  u_num :          numerical solution at grid points 
%  mass, energy :   mass and energy of the numerical solution over time 
%  u_ref :          reference solution    

function [ x, u_num, mass, energy, u_exact, x_eval ] = solve_linear_SAT_inflow_2( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval )

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
    dt = 0.01*(dx_min*block_width); % time-step size
    
    % Global grid points 
    x = zeros(N,I); 
    for i=1:I 
        x(:,i) = x_L + (x_R-x_L)*(i-1)/I + x_ref*block_width;
    end

    % initial data and reference solution 
    u = u_init(x); % solution values on the global grid 
    if x_eval==0 
        x_eval = x; 
    end
    u_exact = u_init(x_eval-T) % reference solution 
   %  u_exact = u_init(x_eval);
    % Weak enforcement of boundary conditions and coupling 
  
    SAT = zeros(N,1);
    sigma = 1;
    mass = []; energy = []; % mass and energy over time 
    
    %% Iterate over time with a 3th-order Runge-Kutta until time T is reached 
    t = 0; 
    while (t<T)  

        % time stepping 
        if T-t<dt 
            dt = T-t; 
        end
        
        t = t+dt 

        % loop over block 
        for i = 1:I

            % Boundary conditions
             % SAT term to weakly enforce the boundary condition
             if i == 1
                g_L = u_init(0.0-t); %Inflow 
            else 
                g_L = u(N,i-1); % Coupling SAT
            end 

            % SAT 
            SAT(1) = -sigma*(u(1,i)-g_L);
            %% FSBP
            % 1st update step 
             % SAT term to weakly enforce the boundary condition
            k1 = u(:,i) + dt*( -D*u(:,i) + P_inv*SAT);  
            % 2nd update step 
          
            k1 = (3/4)*u(:,i) + (1/4)*k1 + (1/4)*dt*( -D*k1 + P_inv*SAT );  
            % 3th update step 
          
            u_num(:,i) = (1/3)*u(:,i) + (2/3)*k1 + (2/3)*dt*( -D*k1 + P_inv*SAT );

        end

        u = u_num; % update solution 
        
        % Compute mass and energy 
        mass_aux = 0; energy_aux = 0; 
        for i=1:I 
            mass_aux = mass_aux + dot( ones(N,1), P*u(:,i) ); % compute mass 
            energy_aux = energy_aux + dot( u(:,i), P*u(:,i) ); % compute energy
        end 
        % save mass and energy
        mass = [mass; t, mass_aux]; % save mass
        energy = [energy; t, energy_aux]; % save energy 
        
    end
    
end

%% Function to compute BCs and SATs 
function [ SAT ] = compute_SAT( u ) 
    
    [ N, I ] = size(u);
    SAT = zeros(N,I); % initialize SAT
    
    for i=1:I 
        % left boundary condition 
        if i==1 
            g_L = u(N,I); 
        else 
            g_L = u(N,i-1); 
        end
        % right boundary condition 
        if i==I 
            g_R = u(1,1); 
        else 
            g_R = u(1,i+1); 
        end
        % SATs 
        SAT(1,i) = -1*( u(1,i) - g_L ); 
        
    end
    
end


