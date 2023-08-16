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

function [ x, u_num, energy, u_exact ] = solve_linear_diffusion_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval,kappa,aa)

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
    dt = 0.01*(dx_min*block_width)^.2; % time-step size
    
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
    
    u_steady = @(x,kappa) ((exp(x/kappa)-1)/(exp(1/(2*kappa))-1));
   % u_exact = u_init(x_eval-T); % reference solution 
     u_exact = u_steady(x_eval, kappa);
    % Weak enforcement of boundary conditions and coupling 
 
    sigma0 = -1; 
    sigma1 = 1;
    energy = []; % mass and energy over time 
    
    g_L=0.0;
    g_R=1.0;
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
          
            % SAT 
            SAT = compute_SAT (u, sigma0, sigma1, kappa, aa, g_L, g_R, D); % SAT term to weakly enforce the boundary condition
            
            %% FSBP
            % 1st update step 
            k1 = u(:,i) + dt*( -aa*D*u(:,i)+ kappa*D*D*u(:,i) + P_inv*SAT(:,i) );  
            % 2nd update step 
            k1 = (3/4)*u(:,i) + (1/4)*k1 + (1/4)* dt*( -aa*D*u(:,i)+ kappa*D*D*u(:,i) + P_inv*SAT );  
            % 3th update step 
            u_num(:,i) = (1/3)*u(:,i) + (2/3)*k1 + (2/3)*dt* dt*( -aa*D*u(:,i)+ kappa*D*D*u(:,i) + P_inv*SAT );

        end

        u = u_num; % update solution 
        
        % Compute energy
        energy_aux = 0; 
        for i=1:I 
  
            energy_aux = energy_aux + dot( u(:,i), P*u(:,i) )+ dot( D*u(:,i), P*D*u(:,i)); % compute energy
        end 
        % save mass and energy
        energy = [energy; t, energy_aux]; % save energy 
        
    end
    
end



%% Function to compute BCs and SATs 
function [ SAT ] = compute_SAT( u, sigma0, sigma1, kappa, aa, g_L, g_R, D) 
    
    [ N, I ] = size(u);
    SAT = zeros(N,I); % initialize SAT
    SAT_L = zeros(N,I);
    SAT_R = zeros(N,I);
    sigma_L = zeros(I);
    sigma_R = zeros(I);
    vx = zeros(N);
    vy = zeros(N);
    vz = zeros(N);
    for i=1:I 
        % left boundary condition 
        vx = D*(u(:,i));
        if i==1 
          SAT_L(1,i)  = sigma0*(u(1,i)-g_L ); 
        else 
            SAT_L(1,i)=0;
         vy =D*(u(:,i-1));
          SAT_L(1,i)  =  sigma_R(1)*(u(0,i)-u(N,i-1))+sigma_R(2)*(vx(1)-vy(N)) ;                          
        end
        % right boundary condition 
        if i==I 
           SAT_R(N,i)  = -sigma1*(u(N,i)-g_R );  
        else
            SAT_R(N,i) =0;
            vz = D*(u(:,i+1));
            SAT_R(N,i) = sigma_L(1)*(u(N,i)-u(0,i+1))+sigma_L(2)*(vx(N)-vz(1));
        end
        % SATs  dfdsf
        SAT(:,i) = SAT_L(:,i)+SAT_R(:,i)
    end
    
end


