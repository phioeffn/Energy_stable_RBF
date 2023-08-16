
%% script_FSBP 
%
% Description: 
% Script to compute an FSBP operator for a given approximation space 
%
% Author: Jan Glaubitz and PHilipp Öffner 
% % Date: Aug 16, 2023
% USed for the spectral convergence investigation 

%% Setting up the script 
clc, clear 

%% Free parameters
x_L = -1.; x_R = 1.; % domain boundaries 
approx_space = 'cubic'; % approximation space (poly, trig, exp, cubic, gauss, multi) 
points = 'equid'; % data points (equid, Halton, random) 
K =  5; % dimension of approximation space
            
%% Basis of F, the corresponding derivatives, a spanning set of G, and moments corresponding to the spanning set of G 
[ basis_F, dx_basis_F, span_G, m_G ] = generate_span( x_L, x_R, approx_space, points, K ); 

%% Compute a positive LS-QF 
if strcmp( approx_space, 'trig')
  	% Compute the trapezoidal rule 
   	N = K+1; 
   	x = linspace(x_L, x_R, N)'; 
   	w = (x_R-x_L)/(N-1)*ones(N,1); 
   	w(1) = 0.5*w(1); w(end) = 0.5*w(end); 
else
  	% Compute the LS rule 
   	[x, w] = compute_LSQF( x_L, x_R, span_G, m_G, points ); % LS-QF
end
    
%% Compute an FSBP operators 
[D, P, Q] = compute_FSBP( basis_F, dx_basis_F, x, w ) % based on LS-QF 
    
