% Compute the objective function
% in : 
%     - m : the estimation of the sequence
%     - P[u/v] : patch approximation of the flow
%     - D[u/v] : the dictionaries
%     - S[u/v] : the sparse approximation of the flow
%     - o  : operators for Fourier transform, TV, Wavelet and optical flow, 
%     - p  : parameters
% out : 
%     - obj : the objective value
%
% Copyright (c) 2019 Timothée Schmoderer


function obj = objective(m, Pu, Pv, Du, Dv, Su, Sv, o, p)
    [Nx, Ny, Nt] = size(m); 
    recon = o.FT  * m - p.data;
    tvxy  = o.TV  * m; 
    xfm   = o.XFM * m; 
    of    = o.OF  * m; 
     
    [u, v]   = unpatch_flow(Pu, Pv, Nx, Ny); 
    [ux, uy] = grad(u); [vx, vy] = grad(v);     

    n_recon = sum(sum(sum(abs(recon).^2))); 
    n_tv    = sum(sum(sum(sqrt(tvxy(:,:,1:2:end).^2 + tvxy(:,:,2:2:end).^2)))); 
    n_xfm   = sum(sum(sum(abs(xfm)))); 
    
    n_of    = sum(sum(sum(abs(of)))); 
    n_sp    = sum(sum(sum(abs(Su)))) + sum(sum(sum(abs(Sv)))); 
    n_tvof  = sum(sum(sum(sqrt(ux.^2 + uy.^2)))) + sum(sum(sum(sqrt(vx.^2 + vy.^2))));  
    
    n_di = 0;
    for t = 1:Nt-1
       n_di   = n_di + sum(sum((Pu(:, :, t) - Du*Su(:, :, t)).^2)) + sum(sum((Pv(:, :, t) - Dv*Sv(:, :, t)).^2));
    end
    
    lambda1 = p.tvWeight; lambda2 = p.xfmWeight; lambda3 = p.ofWeight;
    lambda4 = p.diWeight; lambda5 = p.spWeight;  lambda6 = p.gofWeight;
    
    obj = 0.5 * n_recon + lambda1 * n_tv + lambda2 * n_xfm + lambda3 * n_of + lambda4 * n_di + lambda5 * n_sp + lambda6 * n_tvof;
end

function [ux, uy] = grad(u)
    ux = u([2:end end], :, :) - u; 
    uy = u(:, [2:end end], :) - u;
end