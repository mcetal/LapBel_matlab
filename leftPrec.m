%% Preconditioner for GMRES
%
% The function 
%   y = leftPrec(x, nbk, E, F, Sinv)
% returns the solution to [I | E; F D] y = x
%
% INPUT:
%   x = rhs
%   nbk = total number of unknowns
%   E, F = block matrices for log sources
%   Sinv = inv(D - F E) - the inverse of the Schur complement
% OUTPUT
%   y = [I | E; F | D]^{-1} x
%

function y = leftPrec(x, nbk, E, F, Sinv)
    y = zeros(size(x));
    r_sigma = x(1: nbk);
    r_a = x(nbk+1: end);
    z_a = Sinv*(r_a - F*r_sigma);
    z_sigma = r_sigma - E*z_a;
    y(1:nbk) = z_sigma;
    y(nbk+1: end) = z_a;
end