function [dth, alph, R, T, N, dsda, diagK, Ck] = ...
                                  island_geometry(island_geo, nbod, Np)
                              
    A = island_geo(1, :);
    B = island_geo(2, :);
    th_k =  island_geo(3, :);
    phi_k = island_geo(4, :);
    dth = 2*pi/Np; alph = 0:dth:2*pi-dth;
    
% 
% find coordinate axes so that Zk (=Ck) is through island centres, 
% Xk and Yk are determined accordingly
    [Ck(1,:), Ck(2,:), Ck(3,:)] = sph2cart(th_k, phi_k, ones(1, nbod));
    [Xk(1,:), Xk(2,:), Xk(3,:)] = sph2cart(th_k, phi_k-pi/2, ones(1, nbod));
    Yk = cross(Ck, Xk);
    
    Rp = zeros(3, Np); dRp = Rp; d2Rp = Rp;
    R = zeros(3, nbod*Np); T = R; N = R;
    dsda = zeros(1, nbod*Np);
    diagK = zeros(1, nbod*Np);
    for kbod = 1:nbod
        rotMatrix = [Xk(:, kbod) Yk(:, kbod) Ck(:, kbod)];
        i1 = Np*(kbod-1) + 1; i2 = Np*kbod;
        Rp(1, :) = A(kbod)*cos(alph);
        Rp(2, :) = B(kbod)*sin(alph);
        Rp(3, :) = sqrt(1 - Rp(1, :).^2 - Rp(2, :).^2);
        dRp(1, :) = -A(kbod)*sin(alph);
        dRp(2, :) = B(kbod)*cos(alph);
        dRp(3, :) = (-Rp(1, :).*dRp(1, :) - Rp(2, :).*dRp(2, :))./Rp(3, :); 
        d2Rp(1, :) = -A(kbod)*cos(alph);
        d2Rp(2, :) = -B(kbod)*sin(alph);
        d2Rp(3, :) = (-Rp(1, :).*d2Rp(1, :) - dRp(1, :).^2 ...
                      -Rp(2, :).*d2Rp(2, :) - dRp(2, :).^2 ...
                      -dRp(3, :).^2)./Rp(3, :);
        R(:, i1: i2) = rotMatrix*Rp; 
        T(:, i1: i2) = rotMatrix*dRp;
        dsda(i1: i2) = sqrt(dot(T(:, i1: i2), T(:, i1: i2)));
        T(:, i1: i2) = bsxfun(@rdivide, T(:, i1: i2), dsda(i1: i2));
        N(:, i1: i2) = cross(T(:, i1: i2), R(:, i1: i2));
        dT_ds = bsxfun(@rdivide, rotMatrix*d2Rp, dsda(i1:i2).^2);
        diagK(i1: i2) = -dot(T(:, i1:i2), cross(dT_ds, R(:, i1:i2))) ...
                         .*dsda(i1:i2)/(4*pi);
    end
 end

