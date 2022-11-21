function y = aconv(f,x)
% aconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = aconv(f,x)
%  Inputs
%    f    filter
%    x    1-d signal
%  Outputs
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with the
%    time-reverse of f.
%
%  See Also
%    iconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo

    % dim(x) = Nx x Ny x Nt
    % dim(f) = 1 x 4
    n = size(x, 2); % get Ny 
	p = length(f);
	if p < n
        xpadded = cat(2, x, x(:, (1:p), :)); 
    else
        error("in aconv.m p < n"); 
	   z = zeros(1,p);
	   for i = 1:p
		   imod = 1 + rem(i-1, n);
		   z(i) = x(imod);
	   end
	   xpadded = [x z];
	end
	fflip   = reverse(f);
	ypadded = filter(fflip, 1, xpadded, [], 2); % filter acts on second dimension
          y = ypadded(:, p:(n+p-1), :);
end

    
% 	n = length(x);
% 	p = length(f);
% 	if p < n
% 	   xpadded = [x x(1:p)];
% 	else
% 	   z = zeros(1,p);
% 	   for i=1:p,
% 		   imod = 1 + rem(i-1,n);
% 		   z(i) = x(imod);
% 	   end
% 	   xpadded = [x z];
% 	end
% 	fflip = reverse(f);
% 	ypadded = filter(fflip,1,xpadded);
% 	y = ypadded(p:(n+p-1));

%
% Copyright (c) 1993. David L. Donoho
%     
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 