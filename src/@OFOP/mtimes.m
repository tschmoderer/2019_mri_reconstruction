% in : 
%     - a : OFOP operator 
%     - b : images sequances
% out : 
%     - a*b or a'*b 

function res = mtimes(a, b)
%     mx  = @(m) (m([2:end end],:,:) - m); % derivation selon ---> x
%     my  = @(m) (m(:,[2:end end],:) - m); % dérivation selon ---> y
%     mt  = @(m) (m(:, :, 2:end) - m(:, :, 1:end-1)); % dérivation selon ---> t 
% 
%     yx = @(dm) cat(1, -dm(1,:,:) , dm(1:end-2,:,:) - dm(2:end-1,:,:) , dm(end-1,:,:));
%     yy = @(dm) cat(2, -dm(:,1,:) , dm(:,1:end-2,:) - dm(:,2:end-1,:) , dm(:,end-1,:));
%     yt = @(dm) cat(3, -dm(:,:,1) , dm(:,:,1:end-2) - dm(:,:,2:end-1) , dm(:,:,end-1));
%     
%     mx  = @(m) (m(3:end,:,:) - m(1:end-2,:,:))/2; % derivation selon ---> x
%     my  = @(m) (m(:,3:end,:) - m(:,1:end-2,:))/2; % dérivation selon ---> y
%     mt  = @(m) (m(:,:,2:end) - m(:,:,1:end-1));   % dérivation selon ---> t 
% 
%     yx = @(dm) cat(1, -dm(1,:,:) , dm(1:end-2,:,:) - dm(2:end-1,:,:) , dm(end-1,:,:));
%     yy = @(dm) cat(2, -dm(:,1,:) , dm(:,1:end-2,:) - dm(:,2:end-1,:) , dm(:,end-1,:));
%     yt = @(dm) cat(3, -dm(:,:,1) , dm(:,:,1:end-2) - dm(:,:,2:end-1) , dm(:,:,end-1));

    if a.adjoint 
        res = yt(b) + yx(a.ofU.*b) + yy(a.ofV.*b); 
    else
        res = mt(b) + a.ofU.*mx(b) + a.ofV.*my(b);
    end
end


function y = mt(b)
    y = b(:, :, 2:end) - b(:, :, 1:end-1); 
end

function y = mx(b)
    [Nx, Ny, Nt] = size(b); 
    y = zeros(Nx, Ny, Nt-1); 
    y(2:Nx-1, :, :) = (b(3:Nx, :, 1:Nt-1) - b(1:Nx-2, :, 1:Nt-1))/2; 
end

function y = my(b)
    [Nx, Ny, Nt] = size(b); 
    y = zeros(Nx, Ny, Nt-1); 
    y(:, 2:Ny-1, :) = (b(:, 3:Ny, 1:Nt-1) - b(:, 1:Ny-2, 1:Nt-1))/2; 
end

function b = yt(y)
    [Nx, Ny, Nt] = size(y); 
    Nt = Nt + 1;
    b = zeros(Nx, Ny, Nt); 
    b(:, :, 2:Nt-1) =  y(:, :, 2:Nt-1) - y(:, :, 1:Nt-2); 
    b(:, :, 1)      =  y(:, :, 1); 
    b(:, :, Nt)     = -y(:, :, Nt-1); 
    
    b = -b; 
end

function b = yx(y)
    [Nx, Ny, Nt] = size(y); 
    Nt = Nt + 1;
    b = zeros(Nx, Ny, Nt); 

    b(3:Nx-2, :, 1:Nt-1)  =  y(2:Nx-3, :, :) - y(4:Nx-1, :, :); 
    b(1:2, :, 1:Nt-1)     = -y(2:3, :, :); 
    b(Nx-1:Nx, :, 1:Nt-1) =  y(Nx-2:Nx-1, :, :); 
    b = b/2; 
end

function b = yy(y)
    [Nx, Ny, Nt] = size(y); 
    Nt = Nt + 1;
    b = zeros(Nx, Ny, Nt); 
    
    b(:, 3:Ny-2, 1:Nt-1)  =  y(:, 2:Ny-3, :) - y(:, 4:Ny-1, :); 
    b(:, 1:2, 1:Nt-1)     = -y(:, 2:3, :); 
    b(:, Ny-1:Ny, 1:Nt-1) =  y(:, Ny-2:Ny-1, :);  
    b = b/2; 
end

