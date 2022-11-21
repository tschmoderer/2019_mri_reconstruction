function [um, vm] = resizeFlow2(um, vm, new_size)
    
    scaling = [new_size(1)/size(um,1), new_size(2)/size(um,2)];
    
    a = linspace(1, size(um,2), new_size(2));
    b = linspace(1, size(um,1), new_size(1));
    
    [xi, yi] = meshgrid(a, b);
        
    um = interp2(um.*scaling(2), xi, yi); um(isnan(um)) = 0;
    vm = interp2(vm.*scaling(1), xi, yi); vm(isnan(vm)) = 0;
end