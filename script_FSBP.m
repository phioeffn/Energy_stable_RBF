%% script_FSBP 
%
% Description: 
% Script to compute an FSBP operator for a given approximation space 
%
% Author: J. Glaubitz, J. Nordström and P.Öffner
% Date: Mar 17, 2022 % Author: Jan Glaubitz 


%% Setting up the script 
clc, clear 

%% Free parameters
x_L = 0.; x_R = 1.0; % domain boundaries 
approx_space = 'trig'; % approximation space (poly, trig, exp, cubic, cubic_NEW, gauss, multi) 
points = 'equid'; % data points (equid, Halton, random) 
K =  3; % dimension of approximation space
            
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
    
