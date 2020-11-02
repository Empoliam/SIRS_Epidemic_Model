function [] = drawArrow(t,xt1,xt2,varargin)
%drawArrow Draws an arrow on the given line in the direction of the provided time series
%   Input:
%   
%   t - Column vector, function handles, function to solve
%   xt1 - Column vector, numeric, initial guess
%   xt2 - Function handle, function returning jacobian of f
%   varargin - Optional inputs
%

%Parse varargin
parser = inputParser;

addRequired(parser,'t',@isnumeric);
addRequired(parser,'xt1',@isnumeric);
addRequired(parser,'xt2',@isnumeric);
addOptional(parser,'relDist',0.05,@isnumeric);      %Relative distance along line to draw arrow
addOptional(parser,'lineSpec','k')                  %Line spec
parse(parser,t,xt1,xt2,varargin{:})

relDist = parser.Results.relDist;
lineSpec = parser.Results.lineSpec;

%Determine direction of time series
dir = 1;
if(t(1)>t(2))
    dir = -1;
end

%Determine position of arrow
arcLength = length(xt1);
plotIndex = floor(arcLength.*relDist);

x1 = xt1(plotIndex);
x2 = xt2(plotIndex);

%Calculate arrow direction
dir = [xt1(plotIndex+dir)-xt1(plotIndex), xt2(plotIndex+dir)-xt2(plotIndex)];
dir = dir./norm(dir);

%Draw arrow
arrow3([x1,x2],[x1,x2]+1e-3.*dir,lineSpec,1,1)

end

