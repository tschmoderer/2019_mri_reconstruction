% Compute the errors between the reconstruction and the reference sequence
% in : 
%     - ideals : reference sequence
%     - recons : reconstructed frames
% out : 
%     - err_ssim  : the ssim between each frame
%     - ssim_map  : a map of teh ssim metrics at each pixels
%     - err_slmse : the slmse between each frame
%     - err_rmse  : the rmse between each frame
%     - err_nmse  : the nmse between each frame
%     - err_pnsr  : the pnsr between each frame
%
% Copyright (c) 2019 Timothée Schmoderer

function [err_ssim, ssim_map, err_slmse, err_rmse, err_nmse, err_psnr] = get_errors(ideals, recons)   
    [Nx, Ny, Nt] = size(ideals); 
    
    err_ssim  = zeros(Nt, 1); 
    ssim_map  = zeros(Nx, Ny, Nt); 
    err_slmse = zeros(Nt, 1); 
    err_rmse  = zeros(Nt, 1); 
    err_nmse  = zeros(Nt, 1); 
    err_psnr  = zeros(Nt, 1); 
    
    for t = 1:Nt
        recons(:, :, t) = (recons(:,:,t) - min(min(recons(:,:,t)))) / (max(max(recons(:,:,t))) - min(min(recons(:,:,t))));
        ideals(:, :, t) = (ideals(:,:,t) - min(min(ideals(:,:,t)))) / (max(max(ideals(:,:,t))) - min(min(ideals(:,:,t))));
        
        %%% SSIM
        [err_ssim(t), ssim_map(:,:,t)] = ssim(real(recons(:,:,t)), real(ideals(:,:,t)));
        
        %%% sLMSE
        err_slmse(t) = 1 - lmse(recons(:,:,t), ideals(:,:,t))/lmse(recons(:,:,t), zeros(size(ideals(:,:,t)))); 
        
        %%% RMSE
        [~, err_rmse(t)] = RMSE_cal(ideals(:,:,t), recons(:,:,t));
        
        %%% NMSE
        err_nmse(t ) = mse(recons(:,:,t), ideals(:,:,t)) / mse(zeros(size(recons(:,:,t))), ideals(:,:,t)); 

        %%% PNSR
        err_psnr(t)  = psnr(recons(:,:,t), ideals(:,:,t)); 
    end
end

function [U, rmse] = RMSE_cal(ideal, recon)  
    x = ideal; 
    x_init = recon;

    I = double(abs(x));
    U = double(abs(x_init));
    
    Error = zeros(1:size(U,3)); 

    for ii = 1:size(U,3)
        alpha = sum(dot(U(:,:,ii),I(:,:,ii)))/(sum(dot(U(:,:,ii),U(:,:,ii))));
        U(:,:,ii) = (alpha)*U(:,:,ii);
        E = sum(sum(abs((I(:,:,ii)-U(:,:,ii))).^2));
        E = E*1/(sum(sum(abs((I(:,:,ii))).^2)));
        Error(ii)=(E);
    end

    rmse = mean((Error)); 
end

function S = lmse (tar, ref)
	%%LMSE  Local mean square error
	window_size  = 20;
	window_shift = 10;
    
	[Nx, Ny] = size(ref);

	S = 0;
	for i = 1:window_shift:Nx-window_size
        ix = i:(i+window_size); 
		for j = 1:window_shift:Ny-window_size
            iy = j:(j+window_size);
            S = S + mse(tar(ix, iy), ref(ix, iy));
		end
    end
end

function s = mse(target, ref)
    [Nx, Ny] = size(target); 
    s = sum(sum((ref - target).^2)) / (Nx * Ny); 
end