function x = IWT2_PO(wc, L, qmf)
% IWT2_PO -- Inverse 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT2_PO(wc,L,qmf)
%  Inputs
%    wc    2-d wavelet transform [n by n array, n dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 2d wavelet transform, with
%    wc = FWT2_PO(x,L,qmf), then x = IWT2_PO(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT2_PO, MakeONFilter
    
	[n, J] = quadlength(wc);
	x = wc; 
    
    nc = 2^(L+1);
    for jscal = L:J-1
        top = (nc/2+1):nc; 
        bot = 1:(nc/2); 
        all = 1:nc;
        
        x(all, 1:nc, :) = permute(UpDyadLo(permute(x(bot, 1:nc, :), [2 1 3]), qmf), [ 2 1 3]) + permute(UpDyadHi(permute(x(top, 1:nc, :), [2 1 3]), qmf), [2 1 3]); 
        x(1:nc, all, :) = UpDyadLo(x(1:nc, bot, :), qmf) + UpDyadHi(x(1:nc, top, :), qmf);

        nc = 2*nc;
    end
        
    
end