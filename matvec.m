%% Fast Matrix-Vector Product
% The function y = matvec(x, nbk, nbod, dth, zeta, dzeta, diagK, E, F, D)
% calculates y = [I + K | E; F | D] x using 2-D FMM 
%
% INPUT:
%   x = density and log source strengths
%   nbk = total number of unknowns
%   dth = mesh spacing
%   zeta = source locations in stereographic plane
%   dzeta = d zeta/ d theta
%   diagK = diagonal term in quadrature
%   E, F, D = block matrices for log sources
%   Nlam = number of grid lines in elevation direction
% OUTPUT
%   y = [I + K | E; F | D] x
%

function y = matvec(x, nbk, dth, zeta, dzeta, diagK, E, F, D)
%
% extract of density and solution sources
    y = zeros(size(x,1), 1);
    sigma = x(1: nbk)';
    A_k = x(nbk+1: end);
%
% Set FMM parameters
    iprec = 5;
    nsource = nbk;
    ifcharge = 0;
    source(1,:) = real(zeta);
    source(2,:) = imag(zeta);
    charge = zeros(1, nbk) + 1i*zeros(1, nbk);
    ifdipole = 1;
    dipstr = dth*sigma.*abs(dzeta)/(2*pi);
    dipvec(1,:) = real(-1i*dzeta)./abs(dzeta);
    dipvec(2,:) = imag(-1i*dzeta)./abs(dzeta);
    %tic
    [U] = lfmm2dpart(iprec, nsource, source, ifcharge, charge, ifdipole, ...
                     dipstr,dipvec);
    %disp(['   Time in FMM = ', num2str(toc)])
    if U.ier ~= 0
        disp(['FMM Error Code = ', num2str(U.ier)]);
        return
    end
    zQ = -dth*sigma.*conj(zeta).*dzeta./(1 + abs(zeta).^2);
    zQ = zQ/(2*pi);
    zQsum = imag(sum(zQ));
    y(1:nbk) = 0.5*sigma - real(U.pot) - zQsum + imag(zQ) ...
                         - dth*diagK.*sigma + (E*A_k)';
    y(nbk+1: end) = F*sigma' + D*A_k;

end