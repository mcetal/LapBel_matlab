% The function 
%   [xgrd, ygrd, zgrd, igrd, numgrd] 
%                         = build_grid(nbod, Np, island_geo, Nphi, Nlam)
% constructs grid points on the sphere and indicates which of these
% lie outside of all of the islands
%
% INPUT:
%   nbod = number of islands
%   Np = number of grid points per island
%   island_geo = geometry data for islands
%   Nphi = number of grid lines in azimuthal direction
%   Nlam = number of grid lines in elevation direction
% OUTPUT
%   xgrd, ygrd, zgrd = grid points on sphere
%   igrd(i,j) = 1 if point in fluid domain, 0 otherwise (inside island)
%   numgrd = total number of grid points in fluid domain
%

function [xgrd, ygrd, zgrd, igrd, numgrd] ...
                             = build_grid(nbod, Np, island_geo, Nphi, Nlam)

    A = island_geo(1, :);
    B = island_geo(2, :);
    th_k =  island_geo(3, :);
    phi_k = island_geo(4, :);
    
    [z1, z2, z3] = sph2cart(th_k, phi_k, 1.0);
    z_axis = [z1; z2; z3];
    [x1, x2, x3] = sph2cart(th_k, phi_k-pi/2, 1.0);
    x_axis = [x1; x2; x3];
    y_axis = cross(z_axis,x_axis);

%
% Construct a sphere & determine which points are inside or outside of the
% islands
    dlam = pi/Nlam;  lam =  dlam*(-Nlam: 1: Nlam)';
    dphi = pi/Nphi;  phi =  dphi*(1: 1: Nphi)' - pi/2 - dphi/2;
    [Glam, Gphi] = meshgrid(lam, phi);
    [xgrd, ygrd, zgrd]  = sph2cart(Glam, Gphi, 1.);
    eps = 16*pi/Np;
    [nth, nphi] = size(Glam);
    
    igrd = zeros(nth, nphi);
    for i = 1:nphi
        for j = 1:nth
            point = repmat([xgrd(j, i); ygrd(j, i); zgrd(j, i)], 1, nbod);
            x = dot(point, x_axis); 
            y = dot(point, y_axis);
            z = dot(point, z_axis);
            rad = sqrt((x./A).^2 + (y./B).^2);
            tst = [rad<1+eps; z>0];
            igrd(j, i) = ~any(all(tst, 1)); 
        end
    end
    numgrd = sum(sum(igrd));
    disp(['Number of grid points = ', num2str(numgrd)])
    disp(' ')
end

