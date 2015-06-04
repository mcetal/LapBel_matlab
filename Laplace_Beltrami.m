%% Laplace Beltrami Solver
% The function 
%    [t_direct, t_fmm] = Laplace_Beltrami(island_geo, Np, Nphi, Nlam)
% Solves the Laplace Beltrami equation on a multi-connected domain on the
% unit sphere, via an integral equation formulation. The solution can be
% obtained directly, or by an FMM-accelerated, GMRES solver.
%
% INPUT:
%   island_geo = geometry data for islands,
%       island(k, :) = A_k, B_k, th_k, phi_k
%   Np = number of grid points per island
%   Nphi = number of grid lines in azimuthal direction
%   Nlam = number of grid lines in elevation direction
% OUTPUT
%   t_direct = time taken by direct solver
%   t_fmm = time taken by FMM-accelerated GMRES solver

function [t_direct, t_fmm] = Laplace_Beltrami(island_geo, Np, Nphi, Nlam)
%
% System size
    [~, nbod] = size(island_geo);
    nbk = nbod * Np;
    disp(['Number of Unknowns = ', num2str(nbk + nbod)])
    disp(' ')
    
%
% Define closed curve on S, parametrized by alpha
   [dth, alph, R, T, N, dsda, diagK, Ck] ...
        = island_geometry(island_geo, nbod, Np);
    
%
% Build grid on sphere
    [xgrd, ygrd, zgrd, igrd, numgrd] ...
        = build_grid(nbod, Np, island_geo, Nphi, Nlam);
    
%
% Plot geometry information
    figure(1)
    surf(xgrd, ygrd, zgrd, igrd)
    a=colormap(gray);
    new(1:64,:)=a(64:-1:1,:);
    colormap(new);
    caxis([0 3])
    shading flat
    hold on
    grid off
    axis equal
    
    for kbod = 1: nbod
        i1 = Np*(kbod-1) + 1; i2 = Np*kbod;
        plot3([R(1, i1:i2), R(1, i1)], [R(2, i1:i2), R(2, i1)],...
              [R(3, i1:i2), R(3, i1)], 'k','LineWidth',2)
        plot3(Ck(1, kbod), Ck(2, kbod), Ck(3, kbod), 'k*')
    end
    %quiver3(R(1, :), R(2, :), R(3, :), T(1, :), T(2, :), T(3, :), 0.75)
    %quiver3(R(1, :), R(2, :), R(3, :), N(1, :), N(2, :), N(3, :), 0.75)
    %quiver3(R(1, :), R(2, :), R(3, :), R(1, :), R(2, :), R(3, :), 0.75)
    
    
%
% Get stereographic projection of points
    zeta_k = cart_to_zeta(Ck(1,:), Ck(2,:), Ck(3,:));
    zeta = cart_to_zeta(R(1,:), R(2,:), R(3,:));
    zeta_grd = cart_to_zeta(xgrd, ygrd, zgrd);
    dzeta = (T(1,:) + 1i*T(2,:) + zeta.*T(3,:)).*dsda./(1 - R(3,:));

    figure(2)
        hold on
    for kbod = 1:nbod
        i1 = Np*(kbod-1) + 1;
        i2 = Np*kbod;
        plot (real(zeta(i1: i2)), imag(zeta(i1: i2)))
    end
    title('Stereographic projection')
    xlabel('Real(\zeta)')
    ylabel('Imag(\zeta)')
    
%
% Construct System:  [K | E ] sigma 
%                    [F | D ] A_k   = f
    figure()
    itest = 2;
    f = exact_solution(nbod, zeta_k, zeta, itest);
    f(nbk+1: nbk+nbod) = zeros(nbod, 1);
    for kbod = 1:nbod
        plot (alph, f((kbod-1)*Np+1 :kbod*Np))
        hold on
    end
    title('Boundary Conditions')
    xlabel('\alpha')
    ylabel('f')

    tic
    [K, E, F, D] = build_system(nbod, Np, nbk, dth, R, N, dsda, diagK, Ck);
    SysMat = [K E; F D];
    disp(['Time to build system matrix = ', num2str(toc)])
    disp(' ')
    
%
% Solve directly
    tic
    soln = SysMat\f;
    t_direct = toc;
    disp(['Time to solve directly = ', num2str(t_direct)])
    disp(' ')
%     sigma = soln(1: nbk)';
%     A_k = soln(nbk+1: end)';
    
    
%
% Solve via GMRES/FMM
    Sinv = inv(D - F*E);  % Inverse of Schur complement for preconditioning
    tic;
    tol_gmres = 1.d-10;
    soln_fmm = gmres(@(x) matvec(x, nbk, dth, zeta, dzeta, ...
                                 diagK, E, F, D), f, nbk, tol_gmres, ...
                     nbk, @(x) leftPrec(x, nbk, E, F, Sinv));
    t_fmm = toc;
    disp(['Total GMRES time = ', num2str(t_fmm)])
    sigma = soln_fmm(1: nbk)';
    A_k = soln_fmm(nbk+1:end)';
%
% Debug preconditioner
%     S = [eye(nbk) E; F D];
%     rhs = S*f;
%     f_calc = LeftPrec(rhs, nbk, E, F, Sinv);
%     err_prec = norm(f_calc-f);
%     disp(['Error in Preconditioner = ', num2str(err_prec)])
%     disp(' ')
    
%
% Plot solution
    figure()
    hold on
    for kbod = 1:nbod
        plot (alph, sigma((kbod-1)*Np + 1: Np*kbod))
    end
    title('Layer Density')
    xlabel('\alpha')
    ylabel('\sigma')
    disp(' ')
    disp(['Log source strengths, A_k = ', num2str(A_k)])
    disp(' ')

%
% check solution at a target point
    y_tar = -0.7; x_tar = 0.1;
    z_tar = [x_tar, y_tar, sqrt(1-x_tar^2-y_tar^2)]';
    zeta_tar = cart_to_zeta(z_tar(1), z_tar(2), z_tar(3));
    figure(1)
        plot3(z_tar(1), z_tar(2), z_tar(3), 'r*')
    figure(2)
        plot(real(zeta_tar), imag(zeta_tar), 'r*')
    
    u_calc = double_layer_eval(dth, Ck, R, N, dsda, sigma, A_k, z_tar);
    u_exact = exact_solution(nbod, zeta_k, zeta_tar, itest);
    disp(['Computed solution at target point = ', num2str(u_calc)])
    disp(['   Exact solution at target point = ', num2str(u_exact)])
    disp(['                            Error = ', ...
           num2str(abs(u_calc - u_exact))])
    disp(' ')
    
%
% calculate solution on grid
    tic
    ugrd = solnGrid_FMM(nbk, nbod, dth, numgrd, zeta, dzeta, ...
                        sigma, A_k, Ck, xgrd, ygrd, zgrd, ...
                        igrd, zeta_grd);
    disp(['Time to compute solution on grid by FMM = ', num2str(toc)])
    disp(' ')
    
    
    u_exact = zeros(size(xgrd));
    u_exact(igrd==1) = exact_solution(nbod, zeta_k, zeta_grd(igrd==1), itest);
    err_grd = abs(ugrd - u_exact);
    err_inf = max(max(err_grd));
    disp(['Max error on grid = ', num2str(err_inf)])
    disp(' ')

%
% Plot solution on grid
    figure()
    surf(xgrd, ygrd, zgrd, ugrd)
    shading flat
%   light
   %lighting gourard
    grid off
    axis equal
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    for kbod = 1: nbod
        i1 = Np*(kbod-1) + 1; i2 = Np*kbod;
        plot3([R(1, i1:i2), R(1, i1)], [R(2, i1:i2), R(2, i1)],...
              [R(3, i1:i2), R(3, i1)], 'k','LineWidth',2)
    end
    title('Solution')
%
% plot log error
    figure(6)
    surf(xgrd, ygrd, zgrd,   log10(err_grd))
    shading flat
    grid off
    axis equal
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    for kbod = 1: nbod
        i1 = Np*(kbod-1) + 1; i2 = Np*kbod;
        plot3([R(1, i1:i2), R(1, i1)], [R(2, i1:i2), R(2, i1)],...
              [R(3, i1:i2), R(3, i1)], 'k','LineWidth',2)
    end
    colorbar
    caxis([-16 floor(log10(err_inf))])
    title('Log_{10} Error')

end

%% Unit Sphere to Stereographic Plane
% The function 
%   zeta = cart_to_zeta(x, y, z)
% projects the point (x, y, z) on the unit sphere to the stereographic 
% plane.
%
% INPUT:
%   x, y, z = coordinates of points on unit sphere
% OUTPUT
%   zeta = point in complex plane
%
function zeta = cart_to_zeta(x,y,z)
    zeta = (x + 1i*y)./(1 - z);
end

%% Stereographic Plane to Unit Sphere 
% The function 
%   [x, y, z] = zeta_to_cart(zeta)
% Finds the (x, y, z) on the unit sphere corresponding to zeta in  the  
% stereographic plane.
%
% INPUT:
%   zeta = point in complex plane
% OUTPUT
%   x, y, z = coordinates of points on unit sphere
%
function [x,y,z] = zeta_to_cart(zeta)
    x = (zeta + conj(zeta))./(1+(abs(zeta)).^2);
    y = -1i*(zeta - conj(zeta))./(1+(abs(zeta)).^2);
    z = (zeta*conj(zeta)-1)./(zeta.*conj(zeta)+1);
end