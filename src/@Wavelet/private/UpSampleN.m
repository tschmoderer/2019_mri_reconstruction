function y = UpSampleN(x, s)
% UpSample -- Upsampling operator
%  Usage
%    u = UpSample(d[,s]) 
%  Inputs
%    d   1-d signal, of length n
%    s   upsampling scale, default = 2
%  Outputs
%    u   1-d signal, of length s*n with zeros
%        interpolating alternate samples
%        u(s*i-1) = d(i), i=1,...,n
%

    if nargin == 1
        s = 2;
    end
    N = size(x); 
    
    if length(N) == 2
        N = [N 1];
    end
    
	n = N(2)*s;
	y = zeros(N(1), n, N(3));
	y(:, 1:s:(n-s+1), :) = x;
end
    