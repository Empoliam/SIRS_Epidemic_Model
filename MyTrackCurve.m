function ylist = MyTrackCurve(userf,y0,ytan,varargin)

parser = inputParser;

addRequired(parser,'userf',@(func) isa(func,'function_handle'));
addRequired(parser,'y0',@isnumeric);
addRequired(parser,'ytan',@isnumeric);
addOptional(parser,'stepsize',1e-2,@isnumeric);
addOptional(parser,'maxSolveIter',10,@isnumeric);
addOptional(parser,'nMax',100,@isnumeric);
addOptional(parser,'jStep',1e-6,@isnumeric);
addOptional(parser,'stop',@(y) false,@(func) isa(func,'function_handle'));
addOptional(parser,'sMax',1e-1,@isnumeric);
addOptional(parser,'sMin',1e-12,@isnumeric);
addOptional(parser,'doPrint',false,@islogical);
addOptional(parser,'correctGuess',false,@islogical);
parse(parser,userf,y0,ytan,varargin{:})

stepsize = parser.Results.stepsize;
maxSolveIter = parser.Results.maxSolveIter;
nMax = parser.Results.nMax;
jStep = parser.Results.jStep;
stop = parser.Results.stop;
sMax = parser.Results.sMax;
sMin = parser.Results.sMin;
doPrint = parser.Results.doPrint;
correctGuess = parser.Results.correctGuess;

ylist = NaN(length(y0),nMax);

i = 1;
yk = y0;

while i <= nMax
    
    if (i == 1 && correctGuess)
        stepsize = 0;
    elseif (i == 2 && correctGuess)
        stepsize = parser.Results.stepsize;
    end

    yp = yk + stepsize.*ytan;
    
    f = @(y) [userf(y); ytan.' * (y-yp)];
    df = @(y) MyJacobian(f,y,jStep);
   
    [ykn,converged,~] = Solve(f,yp,df,'tol',1e-6,'maxIter',maxSolveIter);
        
    while(converged == 0)
        
        stepsize = stepsize/2;
        
        if(stepsize < sMin)
            warning("Failed to converge - Step size too small. This issue may potentially be solved by setting correctGuess to false.")
            return
        end
        
        yp = yk + stepsize.*ytan;
        
        f = @(y) [userf(y); ytan.' * (y-yp)];
        df = @(y) MyJacobian(f,y,jStep);
        
        [ykn,converged,~] = Solve(f,yp,df,'tol',1e-6,'maxIter',maxSolveIter);
        
    end
    
    yk = ykn;
    
    if(stop(yk))
        break
    end
    
    ylist(:,i) = yk;
    
     if(stepsize < sMin && i ~= 1)
        break
    end
      
    A = [MyJacobian(userf,yk,jStep);ytan.'];
    I = eye(length(A));
    B = I(:,end);
    
    z=A\B;
    ytan = (z./norm(z,Inf)) * sign(z' * ytan);
                 
    stepsize = min(stepsize*1.25,sMax);
    i = i + 1;
end

if (i < nMax && doPrint)
    disp(strcat("MyTrackCurve converged after ", num2str(i), " iterations. Consider using a smaller value of nMax to save memory."))
end

end

