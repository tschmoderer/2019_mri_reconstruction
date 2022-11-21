function res = mtimes(a, b)
    % res = mtimes(FT, x)
    % res = FT.*x (FT is adjoint or not)
    
    N = size(b); 

    if a.adjoint %% if adjoint computation == inverse fourier transform
        % b is in k-space, res is in image space 
        res = b.*a.mask; 
        res = sqrt(N(1)*N(2))*ifft2(circshift(res, [N(1)/2 N(2)/2 0]));
        res = res.*conj(a.ph); % correct phase       
        
        switch a.mode
            case 0
                res = real(res);
            case 1
                res = real(res);
        end
    else %% normal computation == fourier transform
        switch a.mode
            case 0
                b = real(b);
            case 1
                b = real(b);
        end

        res = b.*a.ph; % phase correct
        res = circshift(fft2(res), [N(1)/2 N(2)/2 0])/sqrt(N(1)*N(2));        
        res = res.*a.mask;
    end

    if size(b,2) == 1
        res = res(:);
    end
end
