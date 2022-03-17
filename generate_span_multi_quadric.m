%% generate_span_cubic
%
% Description: 
%  Generates a basis of F as well as a spanning set of G=(FF)', 
%  where F is an Gaussian RBF approximation space, and the corresponding moments 
% 
% Author: Jan Glaubitz, J. Nordström, and P.Öffner
% Date: Mar. 16, 2022 
%
% INPUT: 
%  x_L, x_R :	domain boundaries 
%  points :     type of center points 
%  K :         	dimension 
%
% OUTPUT: 
%  basis_F :        vector-valued function with basis of F 
%  dx_basis_F :     vector-valued function with derivatives of the basis of F 
%  span_G :         vector-valued function with spanning elements of G
%  m_G :            moments corrsponding to G  

function [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_multi_quadric( x_L, x_R, points, K ) 
 
    % Generate the centers of the RBF interpolants 
    centers = generate_points( points, x_L, x_R, K );
    rbf = @(r) sqrt(1+(r.^2)); % cubic RBF       
    DM = abs(centers-centers');
    A = [ rbf(DM) ones(K,1); ones(1,K) 0 ];
    
    % Coefficients alpha_j, beta_j 
    alpha = zeros(K,K); 
    beta = zeros(1,K); 
    gamma = zeros(K+1,K);
    y = zeros(K+1,K); 
    y(1:K,1:K) = eye(K); 
    gamma = A\y; 
    alpha = gamma(1:K,:); 
    beta = gamma(K+1,:);
       
 	%% Basis of cardinal functions c_k
  	basis_F = @(x) (alpha')*rbf(abs(x-centers)) + beta';  
    syms x 
    dx_basis_F = matlabFunction( diff( basis_F(x) , x ) ); 
    
    %% spanning set for G and corresponding moments 
   	span_G = @(x) basis_F(x)*(dx_basis_F(x)') + dx_basis_F(x)*(basis_F(x)'); 
   	m_G = basis_F(x_R)*(basis_F(x_R)') - basis_F(x_L)*(basis_F(x_L)');
    
  	%% vectorize 
   	span_G = @(x) reshape( span_G(x), [K^2, 1]);
   	m_G = reshape( m_G, [K^2, 1]);
        
end