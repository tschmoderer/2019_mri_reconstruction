function res = adjD(y)
    res = adjDx(y(:,:,(1:2:end))) + adjDy(y(:,:,(2:2:end))); 
end

function res = adjDy(x)
    res            = x(:, [1, 1:end-1], :) - x;
    res(:, 1, :)   = -x(:, 1, :);
    res(:, end, :) = x(:, end-1, :);
end

function res = adjDx(x)
    res            = x([1,1:end-1], :, :) - x;
    res(1, :, :)   = -x(1, :, :);
    res(end, :, :) = x(end-1, :, :);
end