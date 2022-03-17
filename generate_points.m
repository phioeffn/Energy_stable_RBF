%% generate_points
%
% Description: 
%  Function to generate the desired grid points
% 
% Author: Jan Glaubitz, J. Nordström, and P.Öffner
% Date: Mar. 16, 2022 
% 
% INPUT: 
%  points :     String, type of data points 
%  x_L, x_R :   domain boundaries 
%  N :          number of points
%
% OUTPUT: 
%  x : 	vector of points
%       

function x = generate_points( points, x_L, x_R, N )

    x = zeros(1,N);
    
    %% equidistant points 
    if strcmp( points, 'equid')
       	x = linspace(0, 1, N); 
  	
    %% Gauss-Lobatto points 
    elseif strcmp( points, 'Lobatto') 
        [x,w,P]=lglnodes(N-1);
        x = (flip(x)+1)/2; % transform to [0,1] 
        x = (x_R - x_L)*x + x_L; % transform to [x_L,x_R] 
        
    %% Halton points 
    elseif strcmp( points, 'Halton') 
        p = haltonset(1); % generate Halton point set 
       	p = scramble(p,'RR2'); % scramble point set 
       	x = net(p,N)'; % generate the first N points
  	
    %% random points
    elseif strcmp( points, 'random') 
      	rng('default'); rng(1,'twister'); % to make the results reproducable 
      	x = rand(N,1)'; % random points 
  	
    %% else 
    else 
       	error('Desired points not yet implemented!')    
    end
    
    %% Map points from [0,1] and include boundary points
    x = sort(x'); % sort in ascending order
  	x = (x_R-x_L)*x + x_L; % transform from [0,1] to [x_L,x_R]
  	x(1) = x_L; x(end) = x_R; % make sure that boundaries are included
    
end