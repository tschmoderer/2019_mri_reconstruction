function [u, v] = adjoint_patch(Pu, Pv, Nx, Ny)
    [ps2, pn, Nt] = size(Pu); 
    
    ps    = sqrt(ps2);
    wsize = floor(ps/2); 
    

    % center of patches
    ix = 1+wsize : wsize : Nx-wsize+1; 
    iy = 1+wsize : wsize : Ny-wsize+1; 
    
    u = zeros(Nx, Ny, Nt); 
    v = zeros(Nx, Ny, Nt); 
    

    n = 0; 
    for i = ix
        ii = i-wsize:i+wsize-1; 
        for j = iy
            jj = j-wsize:j+wsize-1; 
            n = n + 1;
            u(ii, jj, :) = u(ii, jj, :) + reshape(Pu(:, n, :), ps, ps, Nt); 
            v(ii, jj, :) = v(ii, jj, :) + reshape(Pv(:, n, :), ps, ps, Nt); 
        end
    end
end