% The function 
%   u = double_layer_eval(dth, Ck, R, N, dsda, sigma, A_k, z_tar)
% computes the solution u at a target point z_tar according to 
%   u = dl(sigma) + sum A_k log source
% DIRECTLY
%
% INPUT:
%   dth = mesh spacing
%   Ck = location of log sources (centres of islands)
%   R = boundary points
%   N = principal normal
%   dsda = incremental arclength, ds = dsda*dth;
%   A_k = log source strengths
%   z_tar = target location
% OUTPUT
%   u = solution at z_tar
%

function u = double_layer_eval(dth, Ck, R, N, dsda, sigma, A_k, z_tar)

   dR = bsxfun(@plus, -R, z_tar);
   factor = -dth*sigma.*dsda./(2*pi*dot(dR, dR));
   integrand = factor.*dot(dR, N);
   u = sum(integrand);
   dCk = bsxfun(@plus, -Ck, z_tar);
   logSources = 0.5*A_k.*log(0.5*dot(dCk, dCk));
   u = u + sum(logSources);
   
end
