% reconstruct flow from patches (assume patch are quare)
% [u/v] is size Nx * Ny * Nt-1
% P is size (2 * patch_size * patch_size) * nb_patch * Nt


function [u, v] = unpatch_flow(Pu, Pv, Nx, Ny)
    %% test 
%     Nt = size(Pu, 2); 
%     u = reshape(Pu, [Nx, Ny, Nt]); 
%     v = reshape(Pv, [Nx, Ny, Nt]); 
%     return

%    persistent normal
    [ps2, pn, Nt] = size(Pu); 
    
    ps    = sqrt(ps2);
    wsize = floor(ps/2); 
    
%     Nx = (sqrt(pn) - 1)*wsize + 2*wsize;
%     Ny = Nx;
      
    % center of patches
    ix = 1+wsize : wsize : Nx-wsize+1; 
    iy = 1+wsize : wsize : Ny-wsize+1; 
    
    u = zeros(Nx, Ny, Nt); 
    v = zeros(Nx, Ny, Nt); 
    
%    if isempty(normal)
        normal = zeros(Nx, Ny);
        for i = ix
            ii = i-wsize:i+wsize-1; 
            for j = iy
                jj = j-wsize:j+wsize-1;
                normal(ii, jj) = normal(ii, jj) + 1;  
            end
        end
%    end
    
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

    for t = 1:Nt
        u(:, :, t) = u(:, :, t) ./ normal; 
        v(:, :, t) = v(:, :, t) ./ normal; 
    end
end