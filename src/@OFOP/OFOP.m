classdef OFOP
    properties (Access = public)
        adjoint; 
        ofU; 
        ofV; 
    end
    methods
        function obj = OFOP(Nx, Ny, Nt)
            obj.adjoint = 0; % boolean is adjoint operor (1) or not (0) 
            obj.ofU = zeros(Nx, Ny, Nt - 1); % optical flow U
            obj.ofV = zeros(Nx, Ny, Nt - 1); % optical flow V
        end
        out = ctranspose(obj);
        out = mtimes(obj, arg1); 
    end
end