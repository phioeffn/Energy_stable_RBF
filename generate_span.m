%% generate_span
%
% Description: 
%  Generates a basis of F as well as a spanning set of G=(FF)', where F is the approximation space, and the corresponding moments 
%
% Author: J. Glaubitz, J. Nordström and P.Öffner
% Date: Mar 17, 2022 
% 
% INPUT: 
%  x_L, x_R :       domain boundaries 
%  approx_space :   approximation space F 
%  points :         type of data points (needed for RBF approximation space) 
%  K :              dimension 
%
% OUTPUT: 
%  basis_F :        vector-valued function with basis of F 
%  dx_basis_F :     vector-valued function with derivatives of the basis of F 
%  span_G :         vector-valued function with spanning elements of G
%  m_G :            moments corrsponding to G  

function [ basis_F, dx_basis_F, span_G, m_G ] = generate_span( x_L, x_R, approx_space, points, K )

    if strcmp( approx_space, 'poly')
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_poly( x_L, x_R, K );
    elseif strcmp( approx_space, 'trig')
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_trig( x_L, x_R, K );
    elseif strcmp( approx_space, 'exp') 
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_exp( x_L, x_R, K );  
    elseif strcmp( approx_space, 'cubic') 
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_cubic( x_L, x_R, points, K );
    elseif strcmp( approx_space, 'cubic_NEW') 
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_cubic_NEW( x_L, x_R, points, K );
    elseif strcmp( approx_space, 'gauss') 
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_gauss( x_L, x_R, points, K );
    elseif strcmp( approx_space, 'multi') 
        [ basis_F, dx_basis_F, span_G, m_G ] = generate_span_multi_quadric( x_L, x_R, points, K );
    else
        error('Desired approximation space not yet implemented!')
    end
    
end