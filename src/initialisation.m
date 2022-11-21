% initialise variables;
% in : 
%     - u_data : undersampled data
%     - mask   : the mask 
%     - (optional) m0     : first initialisation of reconstruction
%     - (optional) dataset : string of the name of the dataset used in {"cine", "eth", "phantom_fb", "phantom_bh"}
%     - (optional) acc : acceleration factor in {4, 6, 8}
% out : 
%     - [Nx, Ny, Nt] : size of the data
%     - [Ps, Pn] : patch size and patch number 
%     - Dn : number of atoms in the dictionary 
%     - m0 : first approximation of the sequence obtain by inverse fourier transform
%     - P[u/v]0 : patch of the first approcimation of the optical flow from TV-L1 method
%     - D[u/v]0 : initial dctionary for u and v
%     - S[u/v]0 : initial sparse resentation of u and v 
%     - op0  : the operators for the reconstruction (TF, XFM, OF, TV)
%     - par0 : all the parameters
%     - stat0 : some metrics at the end of initialisation
%
% Copyright (c) 2019 Timothée Schmoderer

function [Nx, Ny, Nt, Ps, Pn, Dn, m0, Pu0, Pv0, Du0, Dv0, Su0, Sv0, op0, par0] = initialisation(u_data, mask, m0, Du, Dv)
    fprintf('Initialize data\n'); 
    tic;
    
    [Nx, Ny, Nt] = size(u_data); 
    Ps = 16; Dn = 1024; 
    Ps = 64; Dn = 64; 

    %% Operators
    op0     = []; 
    op0.FT  = p2DFT(mask, [Nx Ny], 1, 1);  % under sampling fourier operator 
    op0.TV  = TVOP();                      % TV operator 
    op0.XFM = Wavelet('Daubechies', 4, 4); % Wavelet operator 
    op0.OF  = OFOP(Nx, Ny, Nt);            % Optical flow operator

    %% Parameters
    par0           = [];
    par0.data      = u_data; % data undersampled

    par0.Intlim    = 100;      % number of iterations
    par0.threshold = 1e-6;   % error threshold in reconstruction
    
    par0.tvWeight  = 0.003;  % lambda1, default = 0.002
    par0.xfmWeight = 0.0001; % lambda2, default = 0.005
    par0.ofWeight  = 0.001;  % lambda3, default = 0.001
    par0.diWeight  = 0.001;  % lambda4, default = 0.001
    par0.spWeight  = 0.001;  % lambda5
    par0.gofWeight = 0.0001; % lambda6      
    
    %%% parameters for reconstruction 
    par0.reco.Intlim    = 1000;
    par0.reco.threshold = 1e-4;
    par0.reco.tau       = 1/20;
    par0.reco.sigma     = 1/20;
    par0.reco.theta     = 1; 
    
    %%% parameters for sparsity 
    par0.spar.Intlim    = 100;
    par0.spar.sparsity  = floor(0.7*Dn);
    par0.spar.threshold = 1e-4; 
    
    %%% parameters for dico 
    par0.dico.Intlim    = 50;
    par0.dico.threshold = 1e-3; 
    
    %%% parameters for optical flow 
    par0.flow.Intlim    = 200;
    par0.flow.threshold = 1e-3; 
    par0.flow.tau       = 0.25; 
    par0.flow.sigma     = 0.5; 
    par0.flow.theta     = 1; 

    %% initialise images sequence
    if nargin < 3      
       m0 = op0.FT' * u_data;
    end
 
    %% initialise optical flow sequence 
    u0 = zeros(Nx, Ny, Nt-1); v0 = zeros(Nx, Ny, Nt-1); 
    if par0.ofWeight > 0
        for t = 1:Nt-1
           [u0(:, :, t), v0(:, :, t)] = tvl12d(m0(:, :, t), m0(:, :, t+1)); 
        end
    end
    [Pu0, Pv0, Pn] = patch_flow(u0, v0, Ps); 
    
    %% initialize dictionary and sparse representation  
    Su0 = rand(Dn, Pn, Nt-1); Sv0 = rand(Dn, Pn, Nt-1); 
    
    if nargin < 4 || size(Du, 2) ~= Dn || size(Du, 1) ~= Ps*Ps
        if Dn < Pn * (Nt - 1)
            Du0 = Pu0(:, randperm(Pn*(Nt-1), Dn)); 
            Dv0 = Pv0(:, randperm(Pn*(Nt-1), Dn)); 
        else 
            idx = mod(1:Dn, Pn*(Nt-1)-1) + 1; 
            Du0 = Pu0(:, idx); Dv0 = Pv0(:, idx); 
        end
     
        if par0.diWeight > 0
            Du0 = Du0 ./ repmat( sqrt(sum(Du0.^2)), [Ps*Ps, 1] );
            Dv0 = Dv0 ./ repmat( sqrt(sum(Dv0.^2)), [Ps*Ps, 1] );    

            niter_learning = 50; 
            k = 0; 

            no0 = 0; 
            for t = 1:Nt-1
                no0 = no0 + norm(Pu0(:,:,t) - Du0*Su0(:,:,t), 'fro')^2 + norm(Pv0(:,:,t) - Dv0*Sv0(:,:,t), 'fro')^2;
            end    

            errdico = zeros(2*niter_learning, 1); sp = zeros(niter_learning, 1); no_sp = zeros(niter_learning, 1);
            while k < niter_learning
                [Su0, Sv0] = sparsity(Pu0, Pv0, Du0, Dv0, Su0, Sv0, par0);

                sp(k + 1) = (sum(Su0(:) == 0) + sum(Sv0(:) == 0))/(2*length(Su0(:)));
                no_sp(k+1) = sum(abs(Su0(:)) + abs(Sv0(:))); 

                no1 = 0; 
                for t = 1:Nt-1
                    no1 = no1 + norm(Pu0(:,:,t) - Du0*Su0(:,:,t), 'fro')^2 + norm(Pv0(:,:,t) - Dv0*Sv0(:,:,t), 'fro')^2;
                end
                errdico(2*k + 1) = no1; 

                [Du0, Dv0] = dictionary(Pu0, Pv0, Du0, Dv0, Su0, Sv0, par0); 

                no1 = 0; 
                for t = 1:Nt-1
                    no1 = no1 + norm(Pu0(:,:,t) - Du0*Su0(:,:,t), 'fro')^2 + norm(Pv0(:,:,t) - Dv0*Sv0(:,:,t), 'fro')^2;
                end
                errdico(2*k + 2) = no1; 

                if abs(no1 - no0) < 1e-3 * no0
                    break; 
                end
                k = k + 1; 
                no0 = no1; 
            end
        end
    else
       Du0 = Du; Dv0 = Dv;
       if par0.diWeight > 0
           [Su0, Sv0] = sparsity(Pu0, Pv0, Du0, Dv0, Su0, Sv0, par0);
       end      
    end
    
    elaps = toc; 
    fprintf('Initialization successful in %.5f s\n\n', elaps); 
end
