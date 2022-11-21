function [AEE, AE] = get_flow_errors(F, GT)
    u   =  F(:, :, 1); v   =  F(:, :, 2); 
    ugt = GT(:, :, 1); vgt = GT(:, :, 2); 
    
    [Nx, Ny] = size(u); 
    
    AEE = sum(sum(sqrt((u - ugt).^2 + (v - vgt).^2))) / (Nx* Ny); 
    
    U   = cat(3, u,   v,   ones(Nx, Ny)); U   = U   ./ (sqrt(u.^2 + v.^2 + 1)); 
    UGT = cat(3, ugt, vgt, ones(Nx, Ny)); UGT = UGT ./ (sqrt(ugt.^2 + vgt.^2 + 1)); 
    
    AE = sum(sum(acos( dot(U, UGT, 3) ))) / (Nx * Ny); 
end