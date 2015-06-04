function ugrd = solnGrid_FMM(nbk, nbod, dth, numgrd, zeta, ...
                             dzeta, sigma, A_k, Ck, xgrd, ...
                             ygrd, zgrd, igrd, zeta_grd)
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
    ifpot = 0; 
    ifgrad = 0; 
    ifhess = 0;
    ntarget = numgrd;
    target(1, :) = real(zeta_grd(igrd==1));
    target(2, :) = imag(zeta_grd(igrd==1));
    ifpottarg = 1;
    ifgradtarg = 0;
    ifhesstarg = 0;
    %tic
    [U] = lfmm2dpart(iprec, nsource, source, ifcharge, charge, ifdipole, ...
                     dipstr, dipvec, ifpot, ifgrad, ifhess, ntarget, ...
                     target, ifpottarg, ifgradtarg, ifhesstarg);
    %disp(['   Time in FMM = ', num2str(toc)])
    if U.ier ~= 0
        disp(['FMM Error Code = ', num2str(U.ier)]);
        return
    end
    zQsum = imag(sum(-dth*sigma.*conj(zeta).*dzeta./(1 + abs(zeta).^2)));
    zQsum = zQsum/(2*pi);
    potential = - real(U.pottarg) - zQsum;
%
% Add on log sources
    R_tar = zeros(3, ntarget);
    R_tar(1, :) = xgrd(igrd==1);
    R_tar(2, :) = ygrd(igrd==1);
    R_tar(3, :) = zgrd(igrd==1);
    E = zeros(ntarget, nbod);
    for kbod = 1: nbod
        dR = bsxfun(@minus, R_tar, Ck(:, kbod));
        E(:, kbod) = 0.5*log(0.5*dot(dR, dR));
    end
    potential = potential + (E*A_k')';
%
% Calculate solution
    ugrd = zeros(size(xgrd));
    ugrd(igrd==1) = potential;

end
