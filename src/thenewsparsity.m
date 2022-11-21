% minimize over A
% lambda_4*||u-T^{-1}DA||_F^2 + lambda_5*||A||_1
% f2(A) + f1(A)


function [Au1, Av1] = thenewsparsity(m, u, v, Du, Dv, Au0, Av0, p)
    % parameters   
    [Nx, Ny, Nt] = size(m);
    Ps = sqrt(size(Du, 1)); 
    
    lambda4 = p.diWeight; lambda5 = p.spWeight;
    
    [Pu, Pv, Pn] = patch_flow(u, v, Ps); 
    Au1 = zeros(size(Au0)); Av1 = zeros(size(Av0)); 
    
%     DTu = zeros(Dn, Pn, Nt-1); DTv = zeros(Dn, Pn, Nt-1);
%     DDu = 2 * lambda4 * diag(Du' * Du); DDv = 2 * lambda4 * diag(Dv' * Dv); 
%     
%     for t = 1:Nt-1
%        DTu(:,:,t) =  2 * lambda4 * Du' * Pu(:,:,t); 
%        DTv(:,:,t) =  2 * lambda4 * Dv' * Pv(:,:,t); 
%        
%        Au1(:,:,t) = (-lambda5 + DTu) .* (Au0(:,:,t) > 0) + (lambda5 + DTu) .* (Au0(:,:,t) < 0); 
%        Av1(:,:,t) = (-lambda5 + DTv) .* (Av0(:,:,t) > 0) + (lambda5 + DTv) .* (Av0(:,:,t) < 0); 
%        
%        Au1(:,:,t) = Au1(:,:,t) ./ DDu; 
%        Av1(:,:,t) = Av1(:,:,t) ./ DDv; 
%     end
%     
%     return
    
    lambda = 1;   
    gamma_u = lambda5 / (lambda4 * norm(Du'*Du, 'fro') + 10); 
    gamma_v = lambda5 / (lambda4 * norm(Dv'*Dv, 'fro') + 10); 
    k = 1; 
    while k < p.spar.Intlim       
       % smooth step       
       Bu1 = zeros(size(Du, 2), size(Pu, 2), Nt-1); Bv1 = zeros(size(Du, 2), size(Pu, 2), Nt-1);   
       for t = 1:Nt-1
          Bu1(:,:,t) = 2*lambda4/lambda5 * Du'* (Du*Au0(:,:,t) - Pu(:,:,t)); 
          Bv1(:,:,t) = 2*lambda4/lambda5 * Dv'* (Dv*Av0(:,:,t) - Pv(:,:,t));
       end
       
       % [Bu1, Bv1] = gradf2(Au0, Av0, u, v, Du, Dv, Nx, Ny, Nt, Ps, Pn, lambda4); 
       Bu1 = Au0 - gamma_u * Bu1; 
       Bv1 = Av0 - gamma_v * Bv1; 
       
       % proximal step  
       Au1 = (Bu1 + gamma_u) .* (Bu1 < -gamma_u) + (Bu1 - gamma_u) .* (Bu1 > gamma_u); 
       Av1 = (Bv1 + gamma_v) .* (Bv1 < -gamma_v) + (Bv1 - gamma_v) .* (Bv1 > gamma_v); 
       
       % [Au1, Av1] = proxf1(Bu1, Bv1, Du, Dv, Nx, Ny, Nt, Ps, Pn, gamma, lambda5, lambda6); 
       Au1 = Au0 + lambda * (Au1 - Au0);
       Av1 = Av0 + lambda * (Av1 - Av0);
              
       if (norm(Au1(:)-Au0(:), 1) + norm(Av1(:)-Av0(:), 1)) < p.spar.threshold * (norm(Au0(:), 1) + norm(Av0(:), 1))
           break
       end
       
       % iter 
       k = k + 1; 
       Au0 = Au1; Av0 = Av1; 
    end   
end

function [gAu0, gAv0] = gradf2(Au0, Av0, u, v, Du, Dv, Nx, Ny, Nt, Ps, Pn, lambda4)
    tmp_au = zeros(size(Du,1), Pn); tmp_av = zeros(size(Du,1), Pn); 
    tmp_tmp_au = zeros(size(Du,2), Pn); tmp_tmp_av = zeros(size(Du,2), Pn); 

    for t = 1:Nt-1
       tmp_au(:,:,t) = Du*Au0(:,:,t);  
       tmp_av(:,:,t) = Dv*Av0(:,:,t);  
    end
    
    [Pu, Pv, ~] = patch_flow(u, v, Ps); 
    
    for t = 1:Nt-1 
       tmp_tmp_au(:,:,t) = Du'*(Pu(:,:,t) - tmp_au(:,:,t));  
       tmp_tmp_av(:,:,t) = Dv'*(Pv(:,:,t) - tmp_av(:,:,t));  
    end  
    
    gAu0 = -2*lambda4*tmp_tmp_au;
    gAv0 = -2*lambda4*tmp_tmp_av;
    
%     [tmp_au, tmp_av] = adj_unp_unp_flow(tmp_au, tmp_av, Nx, Ny);
%     [tmp_u, tmp_v] = adjoint_unpatch(u, v, Ps);
%     
%     for t = 1:Nt-1 
%        tmp_tmp_au(:,:,t) = Du'*(tmp_u(:,:,t) - tmp_au(:,:,t));  
%        tmp_tmp_av(:,:,t) = Dv'*(tmp_v(:,:,t) - tmp_av(:,:,t));  
%     end
%     
%     gAu0 = -2 * lambda4 * (tmp_tmp_au); 
%     gAv0 = -2 * lambda4 * (tmp_tmp_av);            
end

% chambolle & pock to solve 
% min_C 0.5*||C - B||^2 + gamma*lambda5*||C||_1 + gamma*lambda6||grad(T^{-1}DC||_1
function [Cu1, Cv1] = proxf1(Bu, Bv, Du, Dv, Nx, Ny, Nt, Ps, Pn, gamma, lambda5, lambda6)
    lambda5 = gamma*lambda5; lambda6 = gamma*lambda6;
    
    max_u = Bu >  lambda5; max_v = Bv >  lambda5;
    min_u = Bu < -lambda5; min_v = Bv < -lambda5;
        
    Cu1 = (Bu - lambda5).*max_u + (Bu + lambda5).*min_u;
    Cv1 = (Bv - lambda5).*max_v + (Bv + lambda5).*min_v;
    
    return
    
    tau = 1/1000; sigma = 1/1000; theta = 1; 
    
    Cu0 = zeros(size(Bu)); Cv0 = zeros(size(Bv)); 
    barCu = zeros(size(Bu)); barCv = zeros(size(Bv)); 
    dCu0 = zeros(size(barCu)); dCv0 = zeros(size(barCv));
    dCux0 = zeros(Nx, Ny, Nt-1); dCuy0 = zeros(Nx, Ny, Nt-1);
    dCvx0 = zeros(Nx, Ny, Nt-1); dCvy0 = zeros(Nx, Ny, Nt-1);

    k = 1; 
    while k < 1000
       % dual step 
       [dCu1, dCv1, dCux1, dCuy1, dCvx1, dCvy1] = dual(dCu0, dCv0, dCux0, dCuy0, dCvx0, dCvy0, barCu, barCv, Du, Dv, sigma, Nx, Ny, Nt, Ps, Pn, lambda5, lambda6);
       
       % primal step 
       [Cu1, Cv1] = primal(dCu1, dCv1, dCux1, dCuy1, dCvx1, dCvy1, Cu0, Cv0, Du, Dv, Bu, Bv, tau, Nt, Ps);
       
       % update
        barCu = Cu1 + theta*(Cu1 - Cu0); barCv = Cv1 + theta*(Cv1 - Cv0); 
       
       if (norm(Cu1(:)-Cu0(:)) + norm(Cv1(:)-Cv0(:))) < 1e-3*(norm(Cu0(:)) + norm(Cv0(:)))
           break
       end
       
       % iter
       k = k +1; 
       Cu0 = Cu1; Cv0 = Cv1; 
       dCu0 = dCu1; dCv0 = dCv1;
       dCux0 = dCux1; dCuy0 = dCuy1; dCvx0 = dCvx1; dCvy0 = dCvy1;       
    end
    if k == 1000
        warning('proxf1, k==1000\n');  
    end
end

function [dCu1, dCv1, dCux1, dCuy1, dCvx1, dCvy1] = dual(dCu0, dCv0, dCux0, dCuy0, dCvx0, dCvy0, barCu, barCv, Du, Dv, sigma, Nx, Ny, Nt, Ps, Pn, lambda5, lambda6)
    [dCu1, dCv1, dCux1, dCuy1, dCvx1, dCvy1] = K(barCu, barCv, Du, Dv, Nx, Ny, Nt, Ps, Pn); 
   
    dCu1  = dCu0 + sigma*dCu1; dCv1 = dCv0 + sigma*dCv1; 
    dCux1 = dCux0 + sigma*dCux1; dCuy1 = dCuy0 + sigma*dCuy1; 
    dCvx1 = dCvx0 + sigma*dCvx1; dCvy1 = dCvy0 + sigma*dCvy1;
    
    dCu1 = dCu1 ./ max(1, abs(dCu1/lambda5)); dCv1 = dCv1 ./ max(1, abs(dCv1/lambda5));
    
    no_au = sqrt(dCux1.^2 + dCuy1.^2); no_av = sqrt(dCvx1.^2 + dCvy1.^2); 
    if lambda6 ~= 0
        dCux1 = dCux1 ./ max(1, no_au/lambda6); dCuy1 = dCuy1 ./ max(1, no_au/lambda6);
        dCvx1 = dCvx1 ./ max(1, no_av/lambda6); dCvy1 = dCvy1 ./ max(1, no_av/lambda6);
    else    
        dCux1 = zeros(size(dCux1)); dCuy1 = zeros(size(dCuy1));
        dCvx1 = zeros(size(dCvx1)); dCvy1 = zeros(size(dCvy1));
    end
end

function [Cu1, Cv1] = primal(dCu1, dCv1, dCux1, dCuy1, dCvx1, dCvy1, Cu0, Cv0, Du, Dv, Bu, Bv, tau, Nt, Ps)
    [Cu1, Cv1] = Kt(dCu1, dCv1, dCux1, dCuy1, dCvx1, dCvy1, Du, Dv, Ps, Nt);
    Cu1 = Cu0 - tau * Cu1; 
    Cv1 = Cv0 - tau * Cv1; 
    
    Cu1 = (Cu1 + tau*Bu)/(1 + tau); 
    Cv1 = (Cv1 + tau*Bv)/(1 + tau);   
end

function [dBu, dBv, dBux, dBuy, dBvx, dBvy] = K(Bu, Bv, Du, Dv, Nx, Ny, Nt, Ps, Pn)
    dBu = Bu; dBv = Bv; 
    au = zeros(size(Du,1), Pn, Nt-1); av = zeros(size(Dv, 1), Pn, Nt-1); 
    for t = 1:Nt-1
        au(:,:,t) = Du*Bu(:,:,t);
        av(:,:,t) = Dv*Bv(:,:,t);
    end
    [au, av] = unpatch_flow(au, av, Nx, Ny); 
    
    [dBux, dBuy] = grad(au); 
    [dBvx, dBvy] = grad(av); 
end 

function [Cu, Cv] = Kt(dCu, dCv, dCux, dCuy, dCvx, dCvy, Du, Dv, Ps, Nt)
    du = div(dCux, dCuy); dv = div(dCvx, dCvy); 
    [ddu, ddv] = adjoint_unpatch(du, dv, Ps); 
    
    Cu = dCu; Cv = dCv; 
    for t = 1:Nt-1
       Cu(:, :, t) =  Cu(:,:,t) - Du'*ddu(:,:,t); 
       Cv(:, :, t) =  Cv(:,:,t) - Dv'*ddv(:,:,t); 
    end
end

function [Ux, Uy] = grad(U)
    Ux = U([2:end end], :, :) - U; 
    Uy = U(:, [2:end end], :) - U;
end

function U = div(Ux, Uy)
    U = cat(1, Ux(1, :, :), Ux(2:end-1, :, :) - Ux(1:end-2, :, :), -Ux(end-1, :, :)); 
    U = U + cat(2, Uy(:, 1, :), Uy(:, 2:end-1, :) - Uy(:, 1:end-2, :), -Uy(:, end-1, :)); 
end