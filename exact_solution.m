% The function 
%   f = exact_solution(nbod, zeta_k, zeta, itest)
% returns a function that satisfies Laplace Beltrami at the point(s)
% zeta according to a flag itest. Used to set boundary conditions and check
% computed solution.
%
% INPUT:
%   nbod = number of islands
%   zeta_k = stereographic projection of geometric centres of islands 
%            on sphere
%   zeta = stereographic projection of points along islands
%   itest = a flag for test case
% OUTPUT
%   f = forms rhs for linear system [I + K] mu = f
%
function f = exact_solution(nbod, zeta_k, zeta, itest)

    if itest == 1
        f = real(zeta + conj(zeta));
    elseif itest == 2
        f = zeros(size(zeta));
        for kbod = 1: nbod
            f = f + real(1./(zeta-zeta_k(kbod)) + ...
                          conj(1./(zeta-zeta_k(kbod))));
        end
    end
    f = f';
end

