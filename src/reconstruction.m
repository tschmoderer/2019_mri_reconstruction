% Perform the minization over m of the following energy 
%       0.5*||Km - f||^2 + lambda1*||m||_TV + lambda2*||Psi m||_1 + lambda3*||m_t + grad(m)*U||_1
% in : 
%     - m0 : first estimation of the sequence
%     - o  : operators for Fourier transform, TV, Wavelet and optical flow, 
%     - p  : parameters
% out : 
%     - m1  : new sequence 
%
% Copyright (c) 2019 Timothée Schmoderer

function m1 = reconstruction(m0, o, p) 
    [Nx, Ny, Nt] = size(p.data); 
    
    % Parameters
    m1    = zeros(Nx, Ny, Nt); 
    bar_m = m0; 
    
    sigma = p.reco.sigma;
    tau   = p.reco.tau; 
    theta = p.reco.theta;
    
    Mu = max(abs(o.OF.ofU(:))); Mv = max(abs(o.OF.ofV(:)));
    if tau * sigma > 1 / (2 + sqrt(2) + sqrt(8) + Mu*sqrt(2) + Mv*sqrt(2))^2
       warning('reconstruction: tau and sigma does not satisfied the chambolle and pock condition, tau and sigma decreased');  
       tau   = 0.99 / (2 + sqrt(2) + sqrt(8) + Mu*sqrt(2) + Mv*sqrt(2));
       sigma = 0.99 / (2 + sqrt(2) + sqrt(8) + Mu*sqrt(2) + Mv*sqrt(2));
    end
    
    persistent y1 y2 y3 y4
    if isempty(y1) || isempty(y2) || isempty(3) || isempty(y4)
        y1 = zeros(Nx, Ny, Nt);   % dual of FT
        y2 = zeros(Nx, Ny, 2*Nt); % dual of TV
        y3 = zeros(Nx, Ny, Nt);   % dual of XFM
        y4 = zeros(Nx, Ny, Nt-1); % dual of OF
    end

    % main loop
    k = 1; 
    while (k < p.reco.Intlim)
        tmp_y1 = y1 + sigma * (o.FT  * bar_m); 
        tmp_y2 = y2 + sigma * (o.TV  * bar_m);
        tmp_y3 = y3 + sigma * (o.XFM * bar_m); 
        tmp_y4 = y4 + sigma * (o.OF  * bar_m); 

        % update dual 
        y1 = (tmp_y1 - p.reco.sigma * p.data)/(p.reco.sigma + 1); 
        
        no = sqrt(tmp_y2(:,:,(1:2:end)).^2 + tmp_y2(:,:,(2:2:end)).^2);  
        y2(:,:,(1:2:end)) = tmp_y2(:,:,(1:2:end))./max(1, no / p.tvWeight);
        y2(:,:,(2:2:end)) = tmp_y2(:,:,(2:2:end))./max(1, no / p.tvWeight);
        
        y3 = tmp_y3 ./ max(1, abs(tmp_y3) / p.xfmWeight); 
        y4 = tmp_y4 ./ max(1, abs(tmp_y4) / p.ofWeight);

        % update primal
        tmp_m0 = o.FT'*y1 + (p.tvWeight > 0)*(o.TV'*y2) + (p.xfmWeight > 0)*(o.XFM'*y3) + (p.ofWeight > 0)*(o.OF'*y4); 
        m1     = m0 - tau * tmp_m0;  
        bar_m  = m1 + theta * (m1 - m0);
        
        if norm(m1(:) - m0(:), 2) < p.reco.threshold*norm(m0(:), 2)
            break; 
        end
        
        % iterate
        m0 = m1; 
        k  = k + 1;
    end
    
    if k == p.reco.Intlim
        warning('reconstruction: max iter reached convergence may not be optimal');
    end
end