function res = D(images)
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005

    S = size(images);
    
    if length(S) == 2
        S = [S 1]; 
    end 
    
    res = zeros([S(1) S(2) 2*S(3)]);       
    res(:, :, (1:2:end)) = images([2:end, end] ,:, :) - images;
    res(:, :, (2:2:end)) = images(:, [2:end, end], :) - images;    
end














