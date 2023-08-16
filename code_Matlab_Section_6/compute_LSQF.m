%% compute_LSQF  
%
% Description: 
%  Function to compute the LS-QF points and weights 
%
% Author: Jan Glaubitz 
% Date: Dec 02, 2021 
% 
% INPUT: 
%  x_L, x_R :   domain boundaries 
%  span :       spanning set of G = (FF)'
%  m :          corresponding moments 
%  points :     type of data points 
%
% OUTPUT: 
%  x :          vector of points 
%  w :          vector of weights

function [ x, w] = compute_LSQF( x_L, x_R, span, m, points )

    L = length(m); % number of basis functions (dimension of G)
    N = 1; %1; % dimension of F 
    
    %% routine to determine a nonnegative LS-CF 
    exactness_error = 1; w_min = -1;     
    tol_exactness = 1e-14; % tolerance for the exactness condition    
    while w_min < 1e-14 || exactness_error > tol_exactness
        
        %% data points and matrix G
        x = generate_points( points, x_L, x_R, N ); % points
        G = zeros(L,N); 
        for n=1:N 
            G(:,n) = span(x(n)); 
        end
        
        %% Compute the LS weights 
        w = lsqminnorm(G,m); % indirect computation using optimization tools 
        w_min = min(w); % their smallest value 
        exactness_error = norm( G*w - m )^2/L; 
        
        N = N+1; % increase the number of data points
    
    end
    
end