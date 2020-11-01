clear

sir_model_680029911;
jacobianAccuracy = 1e-6;

load('ylist_stab.mat');
load('hopf_fold.mat');

%Flow map steps
N = 1000;

%Flow map
M = @(t1,x0,p) MyIVP(@(t,x) rhs(x,p),x0,[0,t1],N,'dp45');

%Select hopf bifurcation with largest b value
if(hopfList(3,1) > hopfList(3,2))
    upperHopf = hopfList(:,1);
else
    upperHopf = hopfList(:,2);
end

%Identify closest fold point to hopf point; used to terminate orbit
%tracking
if(foldList(3,1) > foldList(3,2))
    betaStop = foldList(3,1);
else
    betaStop = foldList(3,2);
end

%%Tracking from upper hopf point

%Jacobian of model function
J = @(x) MyJacobian(@(y) rhs(y(1:2),x(3)),x(1:2),jacobianAccuracy);

%Frequency estimate
eVals = eigs(J(upperHopf));
fIni = abs(imag(eVals(1)));

xStar = upperHopf(1);

%[T,x2,p]
xIni = [2.*pi./fIni;upperHopf(2)+1e-2;upperHopf(3)];

%Track orbits until fold point
f = @(x) M(x(1),[xStar;x(2)],x(3)) - [xStar;x(2)];
orbList = MyTrackCurve(f,xIni,[1;0;0],'sMax',0.1,'stop',@(x) x(3) < betaStop,'nMax',150);

%Jacobian of flow map
JM = @(x) MyJacobian(@(y) M(x(1),y,x(3)), [xStar;x(2)], jacobianAccuracy);

%Calculate stability at each orbit, and max/min values
orbStab = NaN(1,length(orbList(1,:)));

%Range of orbits, [beta; Min; Max]
rangeI = NaN(3,length(orbList(1,:)));

i = 1;
while i <= length(orbList(1,:)) && ~isnan(orbList(1,i))
    
    %Calculate stability
    %1 stable,2 unstable
    
    eVals = eigs(JM(orbList(:,i)));
    
    %Handling floating point issues.
    %A tolerance of +-1e-5 is used to identify the eigenvalue equal to 1,
    %such that the remaining eigenvalue can be used to identify stability
    
    if(abs(eVals(1) - 1) < 1e-5)
        compVal =  eVals(2);
    else
        compVal =  eVals(1);
    end
    
    if (compVal <= 1)
        orbStab(i) = 1;
    else
        orbStab(i) = 2;
    end
    
    rangeI(1,i) = orbList(3,i);
    %Calculate min and max    
    %Solve for one full orbit
    [~,~,orb] = MyIVP(@(t,x) rhs(x,orbList(3,i)),[xStar,orbList(2,i)],[0,orbList(1,i)],N,'dp45');
    rangeI(2,i) = min(orb(1,:));
    rangeI(3,i) = max(orb(1,:));
    
    i = i + 1;
end

%Plot bifurcation diagram

cMap = colormap(0.9.*[0,0,1;0,1,0;1,0,0]);

xlabel("beta")
ylabel("I")
xlim([3.5,6]);

plot(ylist(3,:),ylist(1,:),'k-');

hold on
scatter(ylist(3,:),ylist(1,:),15,stab,'filled')
gscatter(rangeI(1,:),rangeI(2,:),orbStab,'gr','__',[5,5],false)
gscatter(rangeI(1,:),rangeI(3,:),orbStab,'gr','++',[5,5],false)
hold off

%Plot orbit family in beta-T plane
figure()

cMap = colormap(0.9.*[0,1,0;1,0,0]);
plot(orbList(3,:),orbList(1,:),'k')
hold on
%Green for stable, red for unstable
scatter(orbList(3,:),orbList(1,:),15,orbStab,'filled')
hold off


%%Plot examples of periodic orbits
figure()

%Select 10 orbits, evenly spaced
orbIndices = floor(linspace(1,nnz(~isnan(orbStab)),10));

i = 1;
while i <= 10
    
    j = orbIndices(i);
    
    %Solve orbit
    [~,~,orb] = MyIVP(@(t,x) rhs(x,orbList(3,j)),[xStar,orbList(2,j)],[0,orbList(1,j)],floor(N/5),'dp45');
    
    %Plot orbit
    subplot(5,2,i)
    
    if(orbStab(j) == 1)
        orbStyle = 'g';
    else
        orbStyle = 'r';
    end
    
    plot(orb(1,:),orb(2,:),orbStyle);
    xlim([0,0.35])
    ylim([0,0.75])
    
    i = i + 1;
    
end