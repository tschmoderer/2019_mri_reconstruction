function y = iconv(f,x)
% iconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv(f,x)
%  Inputs
%    f   filter
%    x   1-d signal
%  Outputs
%    y   filtered result
%
%  Description
%    Filtering by periodic convolution of x with f
%
%  See Also
%    aconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%
	n = size(x, 2); % dim(x) = 1 x Nx x Nt
	p = length(f);
	if p <= n
        xpadded = cat(2, x(:,(n+1-p):n, :), x);
	  % xpadded = [x((n+1-p):n) x];
    else
       error("iconv : p > n"); 
	   z = zeros(1,p);
	   for i=1:p
		   imod = 1 + rem(p*n -p + i-1,n);
		   z(i) = x(imod);
	   end
	   xpadded = [z x];
	end
	ypadded = filter(f, 1, xpadded, [], 2);
          y = ypadded(:, (p+1):(n+p), :);
end

%
% Copyright (c) 1993. David L. Donoho
%     
    
    

    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 