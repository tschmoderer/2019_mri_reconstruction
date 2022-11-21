clc
close all
clear all

%%% 
% - Data are real, 
% - Operators are real
%%%

addpath(strcat(pwd, '/utils')); 
addpath(strcat(pwd, '/utils/flow')); 
addpath(strcat(pwd, '/utils/toolbox_peyre')); 

load data/general_data.mat
ok_save = true; 
verbose = true;

acceleration = acceleration(2:end);

for d = datasets
    load(strcat("data/", d, "_all_cardiac.mat"));
    folder1 = strcat("../results/", d);
    
    for acc = acceleration
        clear reconstruction
        
        load(strcat("../results/", d, "/acceleration_", num2str(acc), "/mcjpdal/", d, "_mcjpdal_", num2str(acc),"_data.mat"), "res"); 
        prev = res;
        
        folder2 = strcat('/acceleration_', num2str(acc), '/dlmcr/', d, '_dlmcr_', num2str(acc));
        
        eval(['mask = '  , char(genvarname(strcat("mask_", num2str(acc)))),';']);
        eval(['u_data = ', char(genvarname(strcat("u_data_", num2str(acc)))),';']);
        
        try
           load(strcat(folder1, folder2, '_data.mat'), 'Du', 'Dv')
        catch
            fprintf("No dictionary found, A new one will be initialized"); 
        end
        
        %% initialisation
        [Nx , Ny, Nt, Ps, Pn, Dn, m0, Pu0, Pv0, Du0, Dv0, Su0, Sv0, op, par] = initialisation(u_data, mask, prev, Du, Dv); 

        %% main function 
        [m1, Pu1, Pv1, Su1, Sv1, stat1] = mri_reconstruction(m0, Pu0, Pv0, Du0, Dv0, Su0, Sv0, op, par, ref); 
        
        [flow_u, flow_v] = unpatch_flow(Pu1, Pv1, Nx, Ny); 
        Du = Du0; Dv = Dv0;
        res = m1;
            
        if verbose
            [ss, ~, ~, ~, ~, ~] = get_errors(ref, res); 
            fprintf("Mean SSIM of %s acceleration factor %d : %.3f\n", d, acc, mean(ss));
        end        
        
        if ok_save
            save(strcat(folder1, folder2, '_data.mat'), 'res', 'flow_u', 'flow_v', 'Du', 'Dv'); 
            for t = 1:Nt
                res(:,:,t) = (res(:,:,t) - min(min(res(:,:,t)))) / (max(max(res(:,:,t))) - min(min(res(:,:,t)))); 
                imwrite(res(:,:,t), strcat(folder1, folder2, '_frame_', num2str(t),'.jpg')); 
            end
            system(strcat("convert -delay 20 -loop 0 ", folder1, folder2, "_frame_*.jpg ", folder1, folder2, '_reconstruction.gif'));
        
            for t = 1:Nt-1
                fl = flowToColor(cat(3, flow_u(:,:,t), flow_v(:,:,t))); 
                imwrite(fl, strcat(folder1,  folder2, '_flow_', num2str(t),'.jpg'));     
            end
            system(strcat("convert -delay 20 -loop 0 ", folder1, folder2, "_flow_*.jpg ", folder1, folder2, '_flow.gif'));
       
            h = imshow_dico(Du0, Dv0); 
            
        end
    end

    if verbose
        fprintf('\n');
    end
    clearvars -except datasets acceleration ok_save verbose
end
