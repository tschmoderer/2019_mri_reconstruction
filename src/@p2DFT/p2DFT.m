function  res = p2DFT(mask, imSize, phase, mode)

%   res = p2DFT(mask, imSize [ , phase, mode])
%
%	Implementation of partial Fourier operator.
%	
%	input:
%			mask   - 2D matrix with 1 in entries to compute the FT and 0 in ones that are not computed.
%			imSize - the image size (1x2)
%			phase  - Phase of the image for phase correction
%			mode   - 0 - positive, 1- real, 2-cmplx
%
%	Output:
%			The operator
%
%	(c) Michael Lustig 2007

    if nargin < 3
        phase = 1;
    end

    if nargin < 4
        mode = 2; % 0 - positive, 1 - real, 2 - cmplx
    end

    res.adjoint  = 0;           % boolean is adjoint operor (1) or not (0) 
    res.mask     = mask;        % the mask 
    res.imSize   = imSize;      % the size of the image 
    res.dataSize = size(mask);  % the size of the fourier data
    res.ph       = phase;       % phase coreection
    res.mode     = mode;        % mode 
    res = class(res,'p2DFT');
end