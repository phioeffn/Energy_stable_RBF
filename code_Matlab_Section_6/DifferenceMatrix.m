% Author: Jan Glaubitz and Philipp Öffner 
% Date: Aug 16, 2021 
% -> Used for the Collocation appraoch

function d = DifferenceMatrix(x)

    [s, N] = size(x);
    d = zeros(N,N);

    for n = 1:N 
        d(n,:) = x(n) - x(:);
    end

end