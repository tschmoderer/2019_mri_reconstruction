clc
clear
close all

addpath(strcat(pwd, "/utils"));
load data/general_data.mat

%% Reference
% test 
%datasets = "cine";
Na = length(acceleration); 
Nd = 6;

zf = 1; cs = 2; cstvl1 = 3; ls = 4; mcjpdal = 5; dlmcr = 6; 

for d = datasets
    SSIM  = zeros(Na, Nd); PSNR  = zeros(Na, Nd); 
    NMSE  = zeros(Na, Nd); RMSE  = zeros(Na, Nd); 
    sLMSE = zeros(Na, Nd); 
    
    load(strcat("data/", d, "_all_cardiac.mat"));
    Nt = size(ref, 3);
    
    for a = 1:Na
        fo = strcat('../results/', d, '/acceleration_', num2str(acceleration(a)));
        
        % ZF
        try
            load(strcat(fo, "/zf/", d, "_zf_", num2str(acceleration(a)), "_data.mat"));
            [ss0, ~, slmse0, rmse0, nmse0, psnr0] = get_errors(ref, res);
        catch
            ss0 = 0; slmse0 = 0; rmse0 = 0; nmse0 = 0; psnr0 = 0;
        end
        
        SSIM(a, zf)  = mean(ss0);   PSNR(a, zf) = mean(psnr0);
        RMSE(a, zf)  = mean(rmse0); NMSE(a, zf) = mean(nmse0);
        sLMSE(a, zf) = mean(slmse0);
        
        % CS
        try
            load(strcat(fo, "/cs/", d, "_cs_", num2str(acceleration(a)), "_data.mat"));
            [ss1, ~, slmse1, rmse1, nmse1, psnr1] = get_errors(ref, res);
        catch
            ss1 = 0; slmse1 = 0; rmse1 = 0; nmse1 = 0; psnr1 = 0;
        end
        
        SSIM(a, cs)  = mean(ss1);   PSNR(a, cs) = mean(psnr1);
        RMSE(a, cs)  = mean(rmse1); NMSE(a, cs) = mean(nmse1);
        sLMSE(a, cs) = mean(slmse1);
        
        % CS+TVL1
        try
            load(strcat(fo, "/cstvl1/", d, "_cstvl1_", num2str(acceleration(a)), "_data.mat"));        
            [ss2, ~, slmse2, rmse2, nmse2, psnr2] = get_errors(ref, res);
        catch 
            ss2 = 0; slmse2 = 0; rmse2 = 0; nmse2 = 0; psnr2 = 0;
        end
        
        SSIM(a, cstvl1)  = mean(ss2);   PSNR(a, cstvl1) = mean(psnr2);
        RMSE(a, cstvl1)  = mean(rmse2); NMSE(a, cstvl1) = mean(nmse2);
        sLMSE(a, cstvl1) = mean(slmse2);
        
        % L+S
        try
            load(strcat(fo, "/ls/", d, "_ls_", num2str(acceleration(a)), "_data.mat"));  
            [ss3, ~, slmse3, rmse3, nmse3, psnr3] = get_errors(ref, res);
        catch
            ss3 = 0; slmse3 = 0; rmse3 = 0; nmse3 = 0; psnr3 = 0;
        end

        SSIM(a, ls)  = mean(ss3);   PSNR(a, ls) = mean(psnr3);
        RMSE(a, ls)  = mean(rmse3); NMSE(a, ls) = mean(nmse3);
        sLMSE(a, ls) = mean(slmse3);        
                
        % MC+JPDAL
        try
            load(strcat(fo, "/mcjpdal/", d, "_mcjpdal_", num2str(acceleration(a)), "_data.mat"));  
            [ss4, ~, slmse4, rmse4, nmse4, psnr4] = get_errors(ref, res);
        catch
            ss4 = 0; slmse4 = 0; rmse4 = 0; nmse4 = 0; psnr4 = 0;
        end

        SSIM(a, mcjpdal)  = mean(ss4);   PSNR(a, mcjpdal) = mean(psnr4);
        RMSE(a, mcjpdal)  = mean(rmse4); NMSE(a, mcjpdal) = mean(nmse4);
        sLMSE(a, mcjpdal) = mean(slmse4);
        
        % DLMC+R
        try
            load(strcat(fo, "/dlmcr/", d, "_dlmcr_", num2str(acceleration(a)), "_data.mat"));  
            [ss5, ~, slmse5, rmse5, nmse5, psnr5] = get_errors(ref, res);
        catch
            ss5 = 0; slmse5 = 0; rmse5 = 0; nmse5 = 0; psnr5 = 0;
        end
        
        SSIM(a, dlmcr)  = mean(ss5);   PSNR(a, dlmcr) = mean(psnr5);
        RMSE(a, dlmcr)  = mean(rmse5); NMSE(a, dlmcr) = mean(nmse5);
        sLMSE(a, dlmcr) = mean(slmse5);
    end
      
    leg  = ["ZF", "CS", "CS+TVL1", "L+S", "MC+JPDAL", "DLMC+R"];
    ttle = strcat("Mean PSNR for the ", d, ' dataset');   
    
    
    
    bar(SSIM)
    % bar(PSNR)
    legend(leg,'FontSize', 18 )

    xticks(1:Na)
    xticklabels({'2x', '4x', '6x', '8x'})

    xlabel('Acceleration factor')
    % ylabel('Mean PSNR')
    
    ylim([0.4 1])
    % ylim([20 40])
    
    title(ttle);
    
    % save(strcat(fo, '/results_', d,'_acc_', acc, '.mat'), 'ss*', 'rmse*', 'psnr*', 'slmse*', 'nmse*');     
    % save(strcat(fo, '/results_', d, '.mat'), 'SSIM', 'PSNR', 'RMSE', 'NMSE', 'sLMSE');

end
