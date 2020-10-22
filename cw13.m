clear

sir_model_680029911;
format long

load('ylist_stab.mat');

f = @(I) rhs(I(1:2),I(3));
df = @(I) MyJacobian(@(a) rhs(a,I(3)),I(1:2),1e-6);

%find approximate bifurcation locations
hopfList = [];
foldList = [];
i = 2;
while i <= length(stab)
    if((stab(i) == 1 && stab(i-1) == 2) || (stab(i) == 2 && stab(i-1) == 1))
        hopfList(end + 1) = i;
    elseif((stab(i) == 0 && stab(i-1) ~= 0) || (stab(i) ~= 0 && stab(i-1) == 0))
        foldList(end + 1) = i;
    end
    i = i + 1;
end

ffold = @(I) [f(I);det(df(I))];
dffold = @(I) MyJacobian(ffold,I,1e-6);

for i = foldList
    
    I0 = ylist(:,i);    
    Solve(ffold,I0,dffold)
    
end

for i = hopfList
    
    I0 = ylist(:,i);
    
    fhopf = @(I) [f(I);trace(df(I))];
    dfhopf = @(I) MyJacobian(fhopf,I,1e-6);
    
    Solve(fhopf,I0,dfhopf)
    
end