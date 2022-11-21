function pdf = createpdf(imSize, p, pctg)
    Nx = max(size(imSize));
    radius = 2; distType = 2; disp = 0;
    while (1)
        try
           pdf = genPDF(imSize, p, pctg, distType, radius, disp);  
           break
        catch
           radius = radius - (1/Nx); 
        end
        if radius < 0
            error('Radius negative, increase p')
        end
    end
end