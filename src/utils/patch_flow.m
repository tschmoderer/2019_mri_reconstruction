% Construct a patch sampling of the flow, 
% Patches are NOT overlapping
% [u/v] is size Nx * Ny * Nt-1
% P[u/v] is size (patch_size * patch_size) * nb_patch * Nt

function [Pu, Pv, nb_patch] = patch_flow(u, v, ps) 
    %% test 
    [Nx, Ny, Nt] = size(u); 
%     Pu = reshape(u, [Nx*Ny, Nt]); 
%     Pv = reshape(v, [Nx*Ny, Nt]); 
%     nb_patch = 1; 
%     return
    
    [Nx, Ny, Nt] = size(u); 
    wsize = floor(ps / 2);
    
    if mod(Nx, ps) ~= 0 || mod(Ny, ps) ~=0
        error("patch_flow.m: patch_size should divide image dimensions"); 
    end
    
    % all the center of patches    
    ix = 1+wsize : wsize : Nx-wsize+1; 
    iy = 1+wsize : wsize : Ny-wsize+1;    
    
    pn = length(ix)*length(iy);     
    Pu = zeros(ps*ps, pn, Nt); 
    Pv = zeros(ps*ps, pn, Nt); 

    nb_patch = 0; 
    for i = ix
        ii = i-wsize:i+wsize-1; 
        for j = iy
            jj = j-wsize:j+wsize-1; 
            nb_patch = nb_patch + 1;  
            Pu(:, nb_patch, :) = reshape(u(ii, jj, :), [ps*ps Nt]); 
            Pv(:, nb_patch, :) = reshape(v(ii, jj, :), [ps*ps Nt]); 
        end
    end
end