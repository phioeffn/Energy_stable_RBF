%% generate_span_trig
%
% Description: 
%  Generates a basis of F as well as a spanning set of G=(FF)', where F is a trigonometric approximation space, and the corresponding moments 
% 
% Author: Jan Glaubitz and Philipp Öffner
% Last change: Dec 02, 2021 
%
% INPUT: 
%  x_L, x_R:	domain boundaries 
%  K:         	dimension 
%
% OUTPUT: 
%  basis_F :        vector-valued function with basis of F 
%  dx_basis_F :     vector-valued function with derivatives of the basis of F 
%  span_G :         vector-valued function with spanning elements of G
%  m_G :            moments corrsponding to G  

function [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_trig( x_L, x_R, K ) 
    
    if rem(K, 2)==0 
        error('Only odd dimensions are allowed!')
    end
        
    %% Especially simple special case
    if K == 1 
        basis_F = @(x) x.^0; 
        dx_basis_F = @(x) 0*x; 
        span_G = @(x) 0*x; 
        m_G = x_R-x_L; 
    
    %% all other cases 
    else  
        %% basis for F and F'
        beta = (1:(K-1)/2)';
        basis_F = @(x) [x.^0; sin( beta*2*pi*(x_R-x_L)*x ); cos( beta*2*pi*(x_R-x_L)*x ) ];
        dx_basis_F = @(x) [0*x; beta*2*pi*(x_R-x_L).*cos( beta*2*pi*(x_R-x_L)*x ); -beta*2*pi*(x_R-x_L).*sin( beta*2*pi*(x_R-x_L)*x ) ];
    
        %% spanning set for G and corresponding moments 
        span_G = @(x) basis_F(x)*(dx_basis_F(x)') + dx_basis_F(x)*(basis_F(x)'); 
        m_G = basis_F(x_R)*(basis_F(x_R)') - basis_F(x_L)*(basis_F(x_L)');
    
        %% vectorize 
        span_G = @(x) reshape( span_G(x), [K^2, 1]);
        m_G = reshape( m_G, [K^2, 1]);
    
    end
    
end