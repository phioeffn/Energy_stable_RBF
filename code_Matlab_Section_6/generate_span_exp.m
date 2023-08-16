%% generate_span_exp
%
% Description: 
%  Generates a basis of F as well as a spanning set of G=(FF)', where F is an exponential approximation space, and the corresponding moments 
% 
% Author: Jan Glaubitz and Philipp Öffner 
% Date: Aug 16, 2023 
%
% INPUT: 
%  x_L, x_R :	domain boundaries 
%  K :         	dimension 
%
% OUTPUT: 
%  basis_F :        vector-valued function with basis of F 
%  dx_basis_F :     vector-valued function with derivatives of the basis of F 
%  span_G :         vector-valued function with spanning elements of G
%  m_G :            moments corrsponding to G  

function [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_exp( x_L, x_R, K ) 
    
    %% Especially simple special cases
    if K == 1 
        error('K needs to be larger than 1!') 
        
    elseif K == 2 
        basis_F = @(x) [ x.^0; exp(x) ]; 
        dx_basis_F = @(x) [ 0*x; exp(x) ]; 
    
    elseif K == 3 
        basis_F = @(x) [ x.^0; x; exp(x) ]; 
        dx_basis_F = @(x) [ 0*x; x.^0; exp(x) ]; 
        
    %% all other cases 
    else       
        %% basis for F and F'
        beta = (0:K-2)'; % exponents 
        basis_F = @(x) [ x.^beta; exp(x) ]; % monomial basis 
        dx_basis_F = @(x) [ 0*x; x.^0; beta(3:end).*x.^(beta(3:end)-1); exp(x) ]; % monomial basis 
    
    end
    
    %% spanning set for G and corresponding moments 
  	span_G = @(x) basis_F(x)*(dx_basis_F(x)') + dx_basis_F(x)*(basis_F(x)'); 
  	m_G = basis_F(x_R)*(basis_F(x_R)') - basis_F(x_L)*(basis_F(x_L)');

  	%% vectorize 
   	span_G = @(x) reshape( span_G(x), [K^2, 1]);
   	m_G = reshape( m_G, [K^2, 1]);
    
end