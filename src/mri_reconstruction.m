% Global function, perform the minimization 
%     0.5*||Km - f||^2 + lambda1*||m||_TV + lambda2*||Psi m||_1 + lambda3*||m_t + grad(m)*U||_1
%                      + lambda4*||TU - DS||^2_F + lambda5*||S||_1 + lambda6*||U||_TV
%     over, m, U, S
% in : 
%     - m0 : first approximation of the sequence                                dim : (Nx)x(Ny)x(Nt)
%     - P[u/v]0 : patch approximation of the optical flow data                  dim : (Ps*Ps)x(Pn)x(Nt-1)
%     - D[u/v]0 : the dictionaries                                              dim : (Ps*Ps)x(Dn)
%     - S[u/v]0 : first approximation of the sparse representation of the flow  dim : (Dn)x(Pn)x(Nt-1)
%     - op  : operators for the reconstruction 
%     - par : parameters
%     - ref : reference images for metrics
% out : 
%     - m1 : reconstructed sequence
%     - P[u/v]1 : new approximation of the flow 
%     - S[u/v]1 : new approximation of the sparse representation of the flow
%     - stat : some convergence metrics
%
% Copyright (c) 2019 Timothée Schmoderer

function [m1, Pu1, Pv1, Su1, Sv1, stat] = mri_reconstruction(m0, Pu0, Pv0, Du, Dv, Su0, Sv0, op, par, ref)
    
    [Nx, Ny, Nt] = size(m0); 
    Ps = sqrt(size(Pu0, 1)); 
    Pn = size(Pu0, 2); 
    Dn = size(Du, 2); 
    
    u1 = zeros(Nx, Ny, Nt-1); v1 = zeros(Nx, Ny, Nt-1); 
    Pu1 = zeros(Ps*Ps, Pn, Nt-1); Pv1 = zeros(Ps*Ps, Pn, Nt-1); 
    Su1 = zeros(Dn, Pn, Nt-1); Sv1 = zeros(Dn, Pn, Nt-1); 
    
    [u0, v0] = unpatch_flow(Pu0, Pv0, Nx, Ny); 

    OBJ = zeros(3*par.Intlim, 1);
	SS  = zeros(par.Intlim, 1);
    for iter = 1:par.Intlim   
        OBJ(3*iter-2) = objective(m0, Pu0, Pv0, Du, Dv, Su0, Sv0, op, par);
        [ss, ~, ~, ~, ~, ~] = get_errors(ref, m0);
        SS(iter) = mean(ss); 
        fprintf('MRI Reconstruction : iteration %03d / %3d  ;  objective : %.3f  ;  mean ssim : %.5f \n', iter, par.Intlim, OBJ(3*iter-2), 100*SS(iter));
        
        if (iter > 2)  && (SS(iter) < SS(iter-1))
            % break
        end
        
        %%%
        % 1 . Reconstruction phase 
        %%%
        op.OF.ofU = u0; op.OF.ofV = v0; 
        m1  = reconstruction(m0, op, par); 
        OBJ(3*iter-2) = objective(m1, Pu0, Pv0, Du, Dv, Su0, Sv0, op, par);   

        %%%
        % 2 . Optical flow phase 
        % 2  . a . Recover optical flow 
        %%%
        if par.ofWeight > 0
             [Su1, Sv1]    = thenewsparsity(m1, u0, v0, Du, Dv, Su0, Sv0, par);
             OBJ(3*iter-1) = objective(m1, Pu0, Pv0, Du, Dv, Su1, Sv1, op, par);

             [u1, v1]      = thenewflow(m1, u0, v0, Du, Dv, Su1, Sv1, par);
             [Pu1, Pv1, ~] = patch_flow(u1, v1, Ps); 
             OBJ(3*iter)   = objective(m1, Pu1, Pv1, Du, Dv, Su1, Sv1, op, par);
        end       

        % iterate 
        m0 = m1; 
        if par.ofWeight > 0, u0 = u1; v0 = v1; Pu0 = Pu1; Pv0 = Pv1; end
        if par.diWeight > 0, Su0 = Su1; Sv0 = Sv1; end
    end

    stat = []; 
    stat.OBJ = OBJ;
    stat.SS  = SS; 
end