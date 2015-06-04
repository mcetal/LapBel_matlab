%% Driver for Laplace_Beltrami.m
%

    close all; clear; clc;
    nbod = 4; 
    A = [0.3   0.2  0.1   0.2];
    B = [-0.2 -0.3 -0.15 -0.25];
    th_k =  [0        3*pi/2   pi/4    pi/2];
    phi_k = [9*pi/20   -pi/4   pi/6  -pi/4];
      
    island_geo = zeros(4, nbod);
    island_geo(1,:) = A(1:nbod);
    island_geo(2,:) = B(1:nbod);
    island_geo(3,:) = th_k(1:nbod);
    island_geo(4,:) = phi_k(1:nbod);
    
    ilev = 6:12;
    t_direct = zeros(1, ilev(end) - ilev(1) + 1);
    t_fmm = zeros(size(t_direct));
    for il = ilev
        close all;  
        Np = 2^il;
        index = il - ilev(1) + 1;
        [t_direct(index), t_fmm(index)] ...
            = Laplace_Beltrami(island_geo, Np, 200, 100);
    end
    
    figure()
    loglog(2.^ilev, t_direct)
    hold on
    loglog(2.^ilev, t_fmm)
    legend('Direct', 'FMM')
    xlabel('N')
    ylabel('Time (seconds)')
    title('Comparison of Direct Solve vs FMM Accelerated')
