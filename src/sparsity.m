% Perform the minization over S of the following energy 
%       ||P0 - D0 * S||_F^2 under contraint ||S_j||_0  <= k0
% in : 
%     - P[u/v]0 : patch of the flow
%     - D[u/v]  : dictionaries reprenstations
%     - S[u/v]0 : first approximation of the sparse reprensetation of the flow
%     - p  : parameters
% out : 
%     - S[u/v]1 : new sparse reprensentation 
%     - err : the convergence error 
%
% Copyright (c) 2019 Timothée Schmoderer

function [Su1, Sv1] = sparsity(Pu0, Pv0, Du, Dv, Su0, Sv0, p)
    Nt1 = size(Su0, 3); 
    Su1 = zeros(size(Su0)); Sv1 = zeros(size(Sv0)); 
    
    gamma_u = 1.6 / norm(Du*Du'); gamma_v = 1.6 / norm(Dv*Dv'); 
    
    for t = 1:Nt1
       [Su1(:, :, t)] = spars(Pu0(:, :, t), Du, Su0(:, :, t), gamma_u, p); 
       [Sv1(:, :, t)] = spars(Pv0(:, :, t), Dv, Sv0(:, :, t), gamma_v, p); 
    end
end

function X1 = ProjX(X0 , k) 
    M = maxk(X0, k); 
    X1 = X0 .* (abs(X0) >= M(end, :));
end

function X1 = spars(P, D, X0, gamma, p) 
    k = 1;
    while k < p.spar.Intlim
        R  = D * X0 - P;
        X1 = ProjX(X0 - gamma * D' * R, p.spar.sparsity);
        
        if sum(abs(X1(:) - X0(:))) < p.spar.threshold * sum(abs(X0(:)))
            break; 
        end
        
        k = k + 1; 
        X0 = X1; 
    end
end