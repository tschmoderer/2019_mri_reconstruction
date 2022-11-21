function res = ctranspose(a)
    % change FT -> FT' and vice versa
    a.adjoint = xor(a.adjoint, 1);
    res = a;
end
