function [] = drawArrow(t,xt1,xt2,varargin)

parser = inputParser;

addRequired(parser,'t',@isnumeric);
addRequired(parser,'xt1',@isnumeric);
addRequired(parser,'xt2',@isnumeric);
addOptional(parser,'relDist',0.05,@isnumeric);
addOptional(parser,'lineSpec','k')
parse(parser,t,xt1,xt2,varargin{:})

relDist = parser.Results.relDist;
lineSpec = parser.Results.lineSpec;

dir = 1;
if(t(1)>t(2))
    dir = -1;
end


arcLength = length(xt1);
plotIndex = floor(arcLength.*relDist);

x1 = xt1(plotIndex);
x2 = xt2(plotIndex);

dir = [xt1(plotIndex+dir)-xt1(plotIndex), xt2(plotIndex+dir)-xt2(plotIndex)];
dir = dir./norm(dir);

arrow3([x1,x2],[x1,x2]+1e-3.*dir,lineSpec,1,1)

end

