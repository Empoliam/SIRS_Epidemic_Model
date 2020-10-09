function [x,converged,J]=Solve(f,x0,df,varargin)

parser = inputParser;

addRequired(parser,'f',@(func) isa(func,'function_handle'));
addRequired(parser,'x0',@isnumeric);
addRequired(parser,'df',@(func) isa(func,'function_handle'));
addOptional(parser,'maxIter',50,@isnumeric);
addOptional(parser,'tol',1e-8,@isnumeric);
parse(parser,f,x0,df,varargin{:})

maxIter = parser.Results.maxIter;
tol = parser.Results.tol;

xn = x0;

iterations = 0;
converged = 0;

while (converged == 0) && (iterations < maxIter)
    
    x = xn;
    J = df(x);
    
    xn = x - (J^(-1)) * f(x);
    
    if (abs(x-xn) < tol)
        converged = 1;
        x = xn;
    end
    
    iterations = iterations + 1;
    
end

end