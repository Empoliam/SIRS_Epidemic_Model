function ylist = MyTrackCurve(userf,userdf,y0,ytan,varargin)

parser = inputParser;

addRequired(parser,'userf',@(func) isa(func,'function_handle'));
addRequired(parser,'userdf',@(func) isa(func,'function_handle'));
addRequired(parser,'y0',@isnumeric);
addRequired(parser,'ytan',@isnumeric);
addOptional(parser,'stepsize',1e-2,@isnumeric);
addOptional(parser,'maxSolveIter',10,@isnumeric);
addOptional(parser,'nMax',100,@isnumeric);
addOptional(parser,'jStep',1e-4,@isnumeric);
addOptional(parser,'stop',@(y) false,@(func) isa(func,'function_handle'));
addOptional(parser,'sMax',1,@isnumeric);
addOptional(parser,'sMin',1e-8,@isnumeric);
addOptional(parser,'doPrint',false,@islogical);
parse(parser,userf,userdf,y0,ytan,varargin{:})

stepsize = parser.Results.stepsize;
maxSolveIter = parser.Results.maxSolveIter;
nMax = parser.Results.nMax;
jStep = parser.Results.jStep;
stop = parser.Results.stop;
sMax = parser.Results.sMax;
sMin = parser.Results.sMin;
doPrint = parser.Results.doPrint;

ylist = zeros(length(y0),nMax);

i = 1;
yk = y0;

while i <= nMax
    
    yp = yk + stepsize.*ytan;
    
    f = @(y) [userf(y); ytan.' * (y-yp)];
    df = @(y) MyJacobian(f,y,jStep);
    
    [ykn,converged,~] = Solve(f,yp,df,'tol',1e-6,'maxIter',maxSolveIter);
    
    while(converged == 0)
        
        stepsize = stepsize/2;
        
        if(stepsize < sMin)
            disp("Failed to converge - Step size too small")
            ykn = NaN(length(yk),1);
            break
        end
        
        yp = yk + stepsize.*ytan;
        
        f = @(y) [userf(y); ytan.' * (y-yp)];
        df = @(y) MyJacobian(f,y,jStep);
        
        [ykn,converged,~] = Solve(f,yp,df,'tol',1e-6,'maxIter',maxSolveIter);
        
    end
    
    yk = ykn;
    
    if(stop(yk))
        i = i -1;
        break
    end
    
    ylist(:,i) = yk;
    
     if(stepsize < sMin)
        break
    end
      
    A = [MyJacobian(userf,yk,jStep);ytan.'];
    I = eye(length(A));
    B = I(:,end);
    
    z=A\B;
    ytan = (z./norm(z,Inf)) * sign(z' * ytan);
     
    if(doPrint)
        disp(i + " - " + "Step Size = " + stepsize )
    end
            
    stepsize = min(stepsize*1.25,sMax);
    i = i + 1;
end

if (i < nMax)
    padding = repmat(ylist(:,i),1,nMax-i+1);
    ylist(:,i:end) = padding;
end

end

