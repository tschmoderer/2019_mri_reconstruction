function [Pu, Pv] = adjoint_unpatch(u, v, ps) 
    %% test 
%     [Nx, Ny, Nt] = size(u); 
%     Pu = reshape(u, Nx*Ny, Nt); 
%     Pv = reshape(v, Nx*Ny, Nt); 
%     return
    
  %  persistent normal
    
    [Nx, Ny, Nt1] = size(u); 
    wsize = floor(ps / 2);
    
    if mod(Nx, ps) ~= 0 || mod(Ny, ps) ~=0
        error("patch_flow.m: patch_size should divide image dimensions"); 
    end
    
    % all the center of patches    
    ix = 1+wsize : wsize : Nx-wsize+1; 
    iy = 1+wsize : wsize : Ny-wsize+1;    
        
%    if isempty(normal)
        normal = zeros(Nx, Ny);
        for i = ix
            ii = i-wsize:i+wsize-1; 
            for j = iy
                jj = j-wsize:j+wsize-1;
                normal(ii, jj) = normal(ii, jj) + 1;  
            end
        end
  %  end

    for t = 1:Nt1
        u(:, :, t) = u(:, :, t) ./ normal; 
        v(:, :, t) = v(:, :, t) ./ normal; 
    end    
    
    pn = length(ix)*length(iy);     
    Pu = zeros(ps*ps, pn, Nt1); 
    Pv = zeros(ps*ps, pn, Nt1); 
    
    nb_patch = 0; 
    for i = ix
        ii = i-wsize:i+wsize-1; 
        for j = iy
            jj = j-wsize:j+wsize-1; 
            nb_patch = nb_patch + 1;  
            Pu(:, nb_patch, :) = reshape(u(ii, jj, :), [ps*ps Nt1]); 
            Pv(:, nb_patch, :) = reshape(v(ii, jj, :), [ps*ps Nt1]); 
        end
    end
end