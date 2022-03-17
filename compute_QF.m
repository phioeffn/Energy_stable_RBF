%% compute_QF  
%
% Description: 
%  Function to compute a positive and (FF)'-exact QF 
%
% Author: Jan Glaubitz, J. Nordström, and P.Öffner
% Date: Mar. 16, 2022 
% 
% INPUT: 
%  x_L, x_R :       domain boundaries 
%  approx_space :   approximation space F 
%  points :         type of data points (needed for RBF approximation space) 
%  K :              dimension 
%
% OUTPUT: 
%  x :	vector of points 
%  w : 	vector of weights

function [ x, w] = compute_QF( x_L, x_R, approx_space, points, K )

    %% Basis of F, the corresponding derivatives, a spanning set of G, and moments corresponding to the spanning set of G 
    [ basis_F, dx_basis_F, span_G, m_G ] = generate_span( x_L, x_R, approx_space, points, K ); 

    %% Compute a positive LS-QF 
    if strcmp( approx_space, 'poly') & strcmp( points, 'Lobatto')
        N = K; 
        [x,w,P]=lglnodes(N-1);
        x = (flip(x)+1)/2; 
        w = w/2;
        
    elseif strcmp( approx_space, 'trig')
        % Compute the trapezoidal rule 
        N = K+1; 
        x = linspace(x_L, x_R, N)'; 
        w = (x_R-x_L)/(N-1)*ones(N,1); 
        w(1) = 0.5*w(1); w(end) = 0.5*w(end); 
    else
        % Compute the LS rule 
        [x, w] = compute_LSQF( x_L, x_R, span_G, m_G, points ); % LS-QF
    end
    
end