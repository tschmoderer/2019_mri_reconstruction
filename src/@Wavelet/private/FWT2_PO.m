function wc = FWT2_PO(x, L, qmf)
% FWT2_PO -- 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT2_PO(x,L,qmf)
%  Inputs
%    x     2-d image (n by n array, n dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    2-d wavelet transform
%
%  Description
%    A two-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT2_PO.
%
%  See Also
%    IWT2_PO, MakeONFilter
	[n, J] = quadlength(x); % n = image size, 2^J = n 
	wc = x; 
    
    nc = n;
    for jscal = J-1:-1:L
        top = (nc/2+1):nc; % index in top row
        bot = 1:(nc/2);    % idex in bottom row 
        
        rows = wc(1:nc, 1:nc, :); 
        wc(1:nc, bot, :) = DownDyadLo(rows, qmf);
        wc(1:nc, top, :) = DownDyadHi(rows, qmf);
        
        rows = permute(wc(1:nc, 1:nc, :), [2 1 3]);
        wc(top, 1:nc, :) = permute(DownDyadHi(rows, qmf), [2 1 3]);
        wc(bot, 1:nc, :) = permute(DownDyadLo(rows, qmf), [2 1 3]); 
        
        nc = nc/2;
    end     
end
