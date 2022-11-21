% Perform the minization over D of the following energy 
%       ||P0 - D * S0||_F^2 under contraint ||D_j||^2  <= 1 
%       perform in parallel minization over the data for u and v 
% in : 
%     - P[u/v]0 : Reference data                     dim : (Ps*Ps)x(Pn)x(Nt)
%     - D[u/v]0 : first estimate of the dictionaries dim : (Ps*Ps)x(Dn)
%     - S[u/v]0 : sparse representation of the flow  dim : (Dn)x(Pn)x(Nt)
%     - p  : parameters
% out : 
%     - D[u/v]1  : new dictionaries
%
% Copyright (c) 2019 Timothée Schmoderer

function [Du1, Dv1] = dictionary(Pu0, Pv0, Du0, Dv0, Su0, Sv0, p)
    Nt1 = size(Pu0, 3);
    
    % preprocess data  
    SSu = 0; PSu = 0;
    SSv = 0; PSv = 0;
    for t = 1:Nt1
       SSu = SSu + Su0(:,:,t)*Su0(:,:,t)'; PSu = PSu + Pu0(:,:,t)*Su0(:,:,t)';
       SSv = SSv + Sv0(:,:,t)*Sv0(:,:,t)'; PSv = PSv + Pv0(:,:,t)*Sv0(:,:,t)';       
    end
    
    tau_u = 1/norm(SSu); tau_v = 1/norm(SSv);
    
    k = 1; 
    while k < p.dico.Intlim
        Ru  = Du0*SSu - PSu;
        Du1 = ProjC(Du0 - tau_u * Ru);
        
        Rv = Dv0*SSv - PSv;
        Dv1 = ProjC(Dv0 - tau_v * Rv);
        
        if norm(Du1 - Du0, 'fro') + norm(Dv1 - Dv0, 'fro') < p.dico.threshold*(norm(Du0, 'fro') + norm(Dv0, 'fro'))
           break;
        end
        
        k = k + 1; 
        Du0 = Du1; Dv0 = Dv1;        
    end
end

function D1 = ProjC(D0)
    D1 = D0 ./ repmat( sqrt(sum(D0.^2)), [size(D0, 1), 1] );
end