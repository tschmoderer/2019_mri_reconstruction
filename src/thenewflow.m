function [u1, v1] = thenewflow(m, u0, v0, Du, Dv, Au, Av, p)
    % parameters   
    [Nx, Ny, Nt] = size(m); Ps = sqrt(size(Du, 1)); Pn = size(Au, 2); 
    lambda3 = p.ofWeight; lambda4 = p.diWeight; lambda6 = p.gofWeight;
    lambda4 = lambda4/lambda3; lambda6 = lambda6/lambda3; lambda3 = 1;
    tau = p.flow.tau; sigma = p.flow.sigma; theta = p.flow.theta;
    
    levels = [4, 2, 1]; 
    fs = fspecial('gaussian', [5 5], 0.8); 
    
    dau = zeros(Ps*Ps, Pn, Nt-1); dav = zeros(Ps*Ps, Pn, Nt-1); 
    for t = 1:Nt-1
       dau(:,:,t) = Du*Au(:,:,t);  
       dav(:,:,t) = Dv*Av(:,:,t); 
    end
    [dau, dav] = adjoint_patch(dau, dav, Nx, Ny); 
    
    % initialisation
    baru = u0; barv = v0; 
    gux0 = zeros(size(u0)); gvx0 = zeros(size(v0)); 
    guy0 = zeros(size(u0)); gvy0 = zeros(size(v0)); 
    
    for gp = levels
       M    = zeros(Nx/gp, Ny/gp, Nt); 
       U0   = zeros(Nx/gp, Ny/gp, Nt-1); V0   = zeros(Nx/gp, Ny/gp, Nt-1);
       GUX0 = zeros(Nx/gp, Ny/gp, Nt-1); GVX0 = zeros(Nx/gp, Ny/gp, Nt-1); 
       GUY0 = zeros(Nx/gp, Ny/gp, Nt-1); GVY0 = zeros(Nx/gp, Ny/gp, Nt-1); 
       BARU = zeros(Nx/gp, Ny/gp, Nt-1); BARV = zeros(Nx/gp, Ny/gp, Nt-1);
       DAU  = zeros(Nx/gp, Ny/gp, Nt-1); DAV  = zeros(Nx/gp, Ny/gp, Nt-1);
       
        wsize = floor((Ps/gp)/2); 
        ix = 1+wsize : wsize : Nx/gp-wsize+1; 
        iy = 1+wsize : wsize : Ny/gp-wsize+1; 

        NORMAL = zeros(Nx/gp, Ny/gp);
        for i = ix
            ii = i-wsize:i+wsize-1; 
            for j = iy
                jj = j-wsize:j+wsize-1;
                NORMAL(ii, jj) = NORMAL(ii, jj) + 1;  
            end
        end

       for t = 1:Nt
          m(:, :, t) = conv2(m(:,:,t), fs, 'same'); 
          M(:, :, t) = imresize(m(:,:,t), [Nx Ny]/gp, 'bilinear'); 
          M(:, :, t) = (M(:,:,t) - min(min(M(:,:,t)))) ./ (max(max(M(:,:,t))) - min(min(M(:,:,t)))); 
          
          if t < Nt
             [U0(:,:,t), V0(:,:,t)]     = resizeFlow2(u0(:,:,t), v0(:,:,t), [Nx Ny]/gp); 
             [BARU(:,:,t), BARV(:,:,t)] = resizeFlow2(baru(:,:,t), barv(:,:,t), [Nx Ny]/gp); 
             [DAU(:,:,t), DAV(:,:,t)]   = resizeFlow2(dau(:,:,t), dav(:,:,t), [Nx Ny]/gp); 
             
             GUX0(:, :, t) = imresize(gux0(:, :, t), [Nx Ny]/gp, 'bilinear'); 
             GUY0(:, :, t) = imresize(guy0(:, :, t), [Nx Ny]/gp, 'bilinear'); 
             GVX0(:, :, t) = imresize(gvx0(:, :, t), [Nx Ny]/gp, 'bilinear'); 
             GVY0(:, :, t) = imresize(gvy0(:, :, t), [Nx Ny]/gp, 'bilinear'); 
          end
       end
       
       [Mt, Mx, My] = dm(M); 
       outer_iter = 1; error = 1; 
       while error > p.flow.threshold && outer_iter < 10
           for t = 1:Nt-1
               U0(:,:,t) = medfilt2(U0(:,:,t), [5 5]); 
               V0(:,:,t) = medfilt2(V0(:,:,t), [5 5]); 
           end
           inner_iter = 1; 
           while error > p.flow.threshold && inner_iter < 30
                % dual step 
                [GUX1, GUY1, GVX1, GVY1] = dual(GUX0, GUY0, GVX0, GVY0, BARU, BARV, sigma, lambda6); 
                
                % primal step 
                [U1, V1] = primal(U0, V0, GUX1, GUY1, GVX1, GVY1, DAU, DAV, Mt, Mx, My, tau, NORMAL, lambda3, lambda4); 
                
                % update 
                BARU = U1 + theta * (U1 - U0); 
                BARV = V1 + theta * (V1 - V0); 
                
                % err
                error = (norm(U1(:)-U0(:)) + norm(V1(:)-V0(:)))/(Nx*Ny*(Nt-1));
                
                %iter 
                inner_iter = inner_iter + 1; 
                U0 = U1; V0 = V1; 
                GUX0 = GUX1; GUY0 = GUY1; 
                GVX0 = GVX1; GVY0 = GVY1; 
           end
           outer_iter = outer_iter + 1; 
       end
         u0 = U1; v0 = V1;
         baru = BARU; barv = BARV; 
         gux0 = GUX1; guy0 = GUY1; 
         gvx0 = GVX1; gvy0 = GVY1; 
    end
    
    u1 = U1; v1 = V1;
end

function [u1, v1] = primal(u0, v0, gux1, guy1, gvx1, gvy1, dau, dav, mt, mx, my, tau, normal, lambda3, lambda4)
	du = div(gux1, guy1); dv = div(gvx1, gvy1); 
    u1 = u0 - tau*du + 2*tau*lambda4*dau; v1 = v0 - tau*dv + 2*tau*lambda4*dav;
    
    A = ones(size(normal)) + 2*tau*lambda4*normal; 
    
    RHO = mt + mx.*(u1 ./ A) + my.*(v1 ./ A); 
    NO  = mx.^2 + my.^2; 
    
    maxus = RHO >  tau*lambda3*NO ./ A; 
    minus = RHO < -tau*lambda3*NO ./ A; 
    other = 1 - maxus - minus;
    
    NO(NO == 0) = 1;
    
    u1 = u1 + (-tau*lambda3*mx).*maxus + (tau*lambda3*mx).*minus + (-RHO.*mx./NO).*other; 
    v1 = v1 + (-tau*lambda3*my).*maxus + (tau*lambda3*my).*minus + (-RHO.*my./NO).*other;    

    u1 = u1 ./ A; 
    v1 = v1 ./ A;     
end

function [gux1, guy1, gvx1, gvy1] = dual(gux0, guy0, gvx0, gvy0, baru, barv, sigma, lambda6)
	[barux, baruy] = grad(baru); 
	[barvx, barvy] = grad(barv);
	
	gux1 = gux0 + sigma*barux; guy1 = guy0 + sigma*baruy;
	gvx1 = gvx0 + sigma*barvx; gvy1 = gvy0 + sigma*barvy;
	
	no_u = sqrt(gux1.^2 + guy1.^2); 
	no_v = sqrt(gvx1.^2 + gvy1.^2); 
	
	gux1 = gux1 ./ max(1, no_u/lambda6); guy1 = guy1 ./ max(1, no_u/lambda6);
	gvx1 = gvx1 ./ max(1, no_v/lambda6); gvy1 = gvy1 ./ max(1, no_v/lambda6);
end

function [mt, mx, my] = dm(m) 
    [Nx, Ny, Nt] = size(m);
    
    mx = zeros(Nx, Ny, Nt - 1); 
    my = zeros(Nx, Ny, Nt - 1); 
    
    mt = m(:,:,2:end) - m(:,:,1:end-1); 
    mx(2:Nx-1, :, :) = (m(3:Nx, :, 1:Nt-1) - m(1:Nx-2, :, 1:Nt-1))/2;
    my(:, 2:Ny-1, :) = (m(:, 3:Ny, 1:Nt-1) - m(:, 1:Ny-2, 1:Nt-1))/2; 
end

function [Ux, Uy] = grad(U)
    Ux = U([2:end end], :, :) - U; 
    Uy = U(:, [2:end end], :) - U;
end

function U = div(Ux, Uy)
    U = cat(1, Ux(1, :, :), Ux(2:end-1, :, :) - Ux(1:end-2, :, :), -Ux(end-1, :, :)); 
    U = U + cat(2, Uy(:, 1, :), Uy(:, 2:end-1, :) - Uy(:, 1:end-2, :), -Uy(:, end-1, :)); 
end