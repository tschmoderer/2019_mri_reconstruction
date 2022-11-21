function d = DownDyadLo(x,qmf)
% DownDyadLo -- Lo-Pass Downsampling operator (periodized)
%  Usage
%    d = DownDyadLo(x,f)
%  Inputs
%    x    1-d signal at fine scale
%    f    filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadHi, UpDyadHi, UpDyadLo, FWT_PO, aconv

    % dim(x) = Nx x Ny x Nt
	d = aconv(qmf, x);
	n = size(d, 2);
	d = d(:, 1:2:(n-1), :);
end