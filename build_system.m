function [K, E, F, D] = build_system(nbod, Np, nbk, dth, R, N, dsda, diagK, Ck)

    K = zeros(nbk);
    E = zeros(nbk, nbod);
    F = zeros(nbod, nbk);
    D = zeros(nbod, nbod);

% Construct discrete integral operator K
    for i = 1: nbk
        Ri = R(:, i);
        dR = bsxfun(@plus, -R, Ri);
        gradU = bsxfun(@rdivide, dR, 2*pi*dot(dR, dR));
        K(i, :) = -dth*dot(gradU, N).*dsda;        
    end
    K(1: nbk+1: nbk^2)  = 0.5-dth*diagK;
    for kbod = 1: nbod
        dR = bsxfun(@minus, R, Ck(:, kbod));
        E(:, kbod) = 0.5*log(0.5*dot(dR, dR));
        i1 = (kbod-1)*Np + 1; i2 = kbod*Np;
        F(kbod, i1: i2) = dsda(i1: i2);
    end
    D(1,:) = ones(1, nbod);
    F(1, 1:Np) = zeros(1, Np);
end