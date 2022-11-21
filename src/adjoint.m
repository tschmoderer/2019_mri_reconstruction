% This script test that operator are well defined and that adjoint are well coded. 

clc
clear
close all

addpath(strcat(pwd, "/utils")); 

% script to test operator and their ajoint
load data/cine_all_cardiac.mat
mask = mask_8;
[Nx, Ny, Nt] = size(data); 

% Operator
FT  = p2DFT(mask, [Nx Ny], 1, 2);  % under sampling fourier operator 
TV  = TVOP();                      % TV operator 
XFM = Wavelet('Daubechies', 4, 4); % Wavelet operator 
OF  = OFOP(Nx, Ny, Nt);            % Optical flow operator

% test operator FT
m = rand(Nx, Ny, Nt); y = rand(Nx, Ny, Nt);
x1 = FT*m; x2 = FT'*y;
fprintf("Adjointness of FT = %e\n", abs(y(:)'*x1(:) - x2(:)'*m(:))); 

% test operator TV 
y = rand(Nx, Ny, 2); m = rand(Nx, Ny);
x1 = TV*m; 
k1 = y(:,:,1)'*x1(:,:,1) + y(:,:,2)'*x1(:,:,2); k2 = (TV'*y)'*m; 
fprintf("Adjointness of TV = %e\n", trace(k1-k2)); 

% test operator XFM
y = rand(Nx, Ny, 2); m = rand(Nx, Ny, 2); 
x1 = XFM*m; x2 = XFM'*y; 
k1 = y(:,:,1)'*x1(:,:,1) + y(:,:,2)'*x1(:,:,2); 
k2 = x2(:,:,1)'*m(:,:,1) + x2(:,:,2)'*m(:,:,2);
fprintf("Adjointness of XFM = %e\n", trace(k1-k2)); 

% test operator OF
Nt = 5; 
m = rand(Nx, Ny, Nt); y = rand(Nx, Ny, Nt-1); 

OF.ofU = rand(Nx,Ny,Nt - 1); 
OF.ofV = rand(Nx,Ny,Nt - 1); 

x1 = OF*m; x2 = OF'*y; 
k1 = 0; k2 = 0; 
for t = 1:Nt-1
    k1 = k1 + y(:,:,t)'*x1(:,:,t); 
end
for t = 1:Nt
    k2 = k2 + x2(:,:,t)'*m(:,:,t);  
end
fprintf("Adjointness of OF = %e\n", trace(k1-k2)); 

% test operator patch 
im1 = rand(Nx, Ny, Nt); 
im2 = rand(Nx, Ny, Nt); 
Ps = gcd(Nx, Ny); 
[Pu, Pv, pn] = patch_flow(im1, im2, Ps); 
[im11, im21] = unpatch_flow(Pu, Pv, Nx, Ny);
fprintf("Invertibility of patch = %e\n", norm(im1(:)-im11(:) + im2(:) - im21(:))); 

yu = rand(size(Pu)); yv = rand(size(Pv)); 
[u1, v1] = adjoint_patch(yu, yv, Nx, Ny); 
fprintf("Adjointness of patch = %e\n", yu(:)'*Pu(:)-u1(:)'*im1(:) + yv(:)'*Pv(:)-v1(:)'*im2(:)); 

% test operator unpatch 
Pu = rand(Ps*Ps, pn, Nt); Pv = rand(Ps*Ps, pn, Nt); 
[u, v] = unpatch_flow(Pu, Pv, Nx, Ny); 
yu = rand(size(u)); yv = rand(size(v)); 
[pu, pv] = adjoint_unpatch(yu, yv, Ps); 

fprintf("Adjointness of unpatch = %e\n", yu(:)'*u(:)-pu(:)'*Pu(:) + yv(:)'*v(:)-pv(:)'*Pv(:)); 

% test operator A in primal of Optical flow 
Ps = 32;
u = rand(Nx, Ny, Nt); v = rand(Nx, Ny, Nt); 

[tu, tv] = patch_flow(u, v, Ps); 
[ttu, ttv] = adjoint_patch(tu, tv, Nx, Ny);

RU1 = u + ttu; RV1 = v + ttv; 

wsize = floor((Ps)/2); 
ix = 1+wsize : wsize : Nx-wsize+1; 
iy = 1+wsize : wsize : Ny-wsize+1; 

normal = zeros(Nx, Ny);
for i = ix
    ii = i-wsize:i+wsize-1; 
    for j = iy
        jj = j-wsize:j+wsize-1;
        normal(ii, jj) = normal(ii, jj) + 1;  
    end
end
A = ones(size(normal)) + normal;

RU2 = A.*u; RV2 = A.*v; 
u1 = RU1./A; v1 = RV1./A; 

fprintf("Definition of A = %e\n", norm(RU1(:) + RV1(:) - RU2(:) - RV2(:))); 
fprintf("Invertibility of A = %e\n", norm(u(:) + v(:) - u1(:) - v1(:))); 

