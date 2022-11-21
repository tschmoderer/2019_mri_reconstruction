function y = UpDyadLo(x,qmf)
% UpDyadLo -- Lo-Pass Upsampling operator; periodized
%  Usage
%    u = UpDyadLo(d,f)
%  Inputs
%    d    1-d signal at coarser scale
%    f    filter
%  Outputs
%    u    1-d signal at finer scale
%
%  See Also
%    DownDyadLo, DownDyadHi, UpDyadHi, IWT_PO, iconv
%
	y =  iconv(qmf, UpSampleN(x));
end 