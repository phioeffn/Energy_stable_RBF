%% compute_FSBP
%
% Description: 
%  Function to compute FSBP operators 
%
% Author: J. Glaubitz, J. Nordström and P.Öffner
% Date: Mar 17, 2022 
% 
% INPUT: 
%  basis_F :        basis of the approximation space F
%  dx_basis_F :     derivatives of the basis elements 
%  x :              grid points 
%  w :              weights of a positive and (F^2)'-exact QF 
%
% OUTPUT: 
%  D :  differentiation matrix 
%  P :  diagonal-norm matrix
%  Q :  matrix for boundary correction         

function [D, P, Q] = compute_FSBP( basis_F, dx_basis_F, x, w )

    N = length(x); % number of grid points 
    K = length( basis_F(x(1)) ); % dimension of F 

    %% Diagonal-norm matrix P 
    P = sparse(diag(w)); 

    %% Anti-symmetric part of Q, Q_anti
    % Prepare matrices 
    F = zeros(N,K); % Vandermonde-like matrix 
    F_x = zeros(N,K); % Vandermonde-like matrix for the derivatives 
    for n=1:N 
    	F(n,:) = basis_F( x(n) )'; 
        F_x(n,:) = dx_basis_F( x(n) )'; 
    end 
    B = zeros(N); B(1,1) = -1; B(end,end) = 1; % boundary matrix  
    R = P*F_x - 0.5*B*F; % right-hand side of the matrix equation for Q_anti 
    
    % Vectorize 
    A = kron( F', eye(N) ); % coefficient matrix for vectorized version 
    r = R(:); % vectorized right-hand side 
    
    % Commutation matrix C
    I = reshape(1:N*N, [N, N]); % initialize a matrix of indices 
    I = I'; % transpose it
    I = I(:); % vectorize the required indices
    C = speye(N*N); % Initialize an identity matrix
    C = C(I,:); % Re-arrange the rows of the identity matrix
    
    % Get the anti-symmetric least-squares solution to 
    A_ext = [ A; C+speye(N^2) ]; 
    r_ext = [ r; zeros(N^2,1) ]; 
    q_anti = lsqminnorm(A_ext,r_ext); % least-squares solution  
    Q_anti = reshape(q_anti,N,N); 
    
    %% Q and the differentiation matrix D 
    Q = Q_anti + 0.5*B; % matrix Q 
    P_inv = diag(1./w); % inverse of P 
    D = P_inv*Q; % differentiation matrix D 
    
end