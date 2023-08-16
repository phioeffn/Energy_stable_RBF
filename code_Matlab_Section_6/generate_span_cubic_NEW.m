%% generate_span_cubic_NEW
%
% Description: 
%  Generates a basis of F as well as a spanning set of G=(FF)', 
%  where F is an cubic RBF approximation space, and the corresponding moments 
% 
% Author: Jan Glaubitz and Philipp Öffner 
% Date: Aug 16, 2023 
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

function [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_cubic_NEW( x_L, x_R, points, K ) 
 
    % Generate the centers of the RBF interpolants 
    centers = generate_points( points, x_L, x_R, K );
    rbf = @(r) r.^3; % cubic RBF       
    basis_F = @(x) [ 1; rbf( abs( x'-centers ) ) ]; % basis of a larger space 
    L = length(basis_F(centers(1))); 
    syms x 
    dx_basis_F = matlabFunction( diff( basis_F(x) , x ) ); 
    
    %% spanning set for G and corresponding moments 
   	span_G = @(x) basis_F(x)*(dx_basis_F(x)') + dx_basis_F(x)*(basis_F(x)'); 
   	m_G = basis_F(x_R)*(basis_F(x_R)') - basis_F(x_L)*(basis_F(x_L)');
    
  	%% vectorize 
   	span_G = @(x) reshape( span_G(x), [L^2, 1]);
   	m_G = reshape( m_G, [L^2, 1]);
        
end