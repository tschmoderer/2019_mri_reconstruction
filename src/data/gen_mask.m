% Generate the undersampling masks
%
% Need:
%     - cine.mat
%     - eth.mat
%     - phantom_bh.mat
%     - phantom_fb.mat
% which can be found in data_mri.tar.xz
% 
% 
% Copyright (c) 2019 Timothée Schmoderer

clc
clear
close all

addpath('./utils/'); 

Nx = 128; Ny = 128; 

imSize   = [Nx 1]; 
nb_try   = 5000; 
nb_err   = 0.5;

% %%% mask 1
% pctg   = 1; p = 1;
% pdf    = createpdf(imSize, p, pctg);
% mask_1 = genSampling(pdf, nb_try, nb_err);
% mask_1 = repmat(mask_1, Ny, 1);
% figure, imshow(mask_1); 
% fprintf('acceleration factor ; %.3f\n',  length(mask_1(:))/ sum(mask_1(:) ~= 0));

%%% mask 2
pctg   = 1/2; p = 3;
pdf    = createpdf(imSize, p, pctg);
mask_2 = genSampling(pdf, nb_try, nb_err);
mask_2 = repmat(mask_2, Ny, 1);
figure, imshow(mask_2); title(["Mask 2 :  ", num2str(pctg*100) , "% of coeffs"]);
fprintf('acceleration factor ; %.3f\n',  length(mask_2(:))/ sum(mask_2(:) ~= 0));

%%% mask 4
pctg   = 1/4; p = 5;
pdf    = createpdf(imSize, p, pctg);
mask_4 = genSampling(pdf, nb_try, nb_err);
mask_4 = repmat(mask_4, Ny, 1);
figure, imshow(mask_4); title(["Mask 4 :  ", num2str(pctg*100) , "% of coeffs"]);
fprintf('acceleration factor ; %.3f\n',  length(mask_4(:))/ sum(mask_4(:) ~= 0));

%%% mask 6
pctg   = 1/6; p = 5;
pdf    = createpdf(imSize, p, pctg);
mask_6 = genSampling(pdf, nb_try, nb_err);
mask_6 = repmat(mask_6, Ny, 1);
figure, imshow(mask_6); title(["Mask 6 :  ", num2str(pctg*100) , "% of coeffs"]);
fprintf('acceleration factor ; %.3f\n',  length(mask_6(:))/ sum(mask_6(:) ~= 0));

%%% mask 8
pctg   = 1/8; p = 7;
pdf    = createpdf(imSize, p, pctg);
mask_8 = genSampling(pdf, nb_try, nb_err);
mask_8 = repmat(mask_8, Ny, 1);
figure, imshow(mask_8); title(["Mask 8 :  ", num2str(pctg*100) , "% of coeffs"]);
fprintf('acceleration factor ; %.3f\n',  length(mask_8(:))/ sum(mask_8(:) ~= 0));

% %%% mask 10
% pctg    = 1/10; p = 9;
% pdf     = createpdf(imSize, p, pctg);
% mask_10 = genSampling(pdf, nb_try, nb_err);
% mask_10 = repmat(mask_10, Ny, 1);
% figure, imshow(mask_10); 
% fprintf('acceleration factor ; %.3f\n',  length(mask_10(:))/ sum(mask_10(:) ~= 0));
% 
% %%% mask 12
% pctg    = 1/12; p = 11;
% pdf     = createpdf(imSize, p, pctg);
% mask_12 = genSampling(pdf, nb_try, nb_err);
% mask_12 = repmat(mask_12, Ny, 1);
% figure, imshow(mask_12); 
% fprintf('acceleration factor ; %.3f\n',  length(mask_12(:))/ sum(mask_12(:) ~= 0));

save mask.mat 'mask_*'

