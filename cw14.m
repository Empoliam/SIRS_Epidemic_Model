clear

sir_model_680029911;
jacobianAccuracy = 1e-8;

load('ylist_stab.mat');
load('hopf_fold.mat');

%Flow map steps
N = 500;

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
fIni = abs(imag(eVals(1))); %Imaginary part of positive eigenvalue

%Fix I at the I value of the hopf point
xStar = upperHopf(1);

%[T,x2,p] initial guess for orbit
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
    
    %Store beta values, for easier plotting
    rangeI(1,i) = orbList(3,i);
    
    %Solve for one full orbit
    [~,~,orb] = MyIVP(@(t,x) rhs(x,orbList(3,i)),[xStar,orbList(2,i)],[0,orbList(1,i)],N,'dp45');
    
    %Calculate min and max
    rangeI(2,i) = min(orb(1,:));
    rangeI(3,i) = max(orb(1,:));
    
    i = i + 1;
end

%Plot bifurcation diagram

figure(1)
cMap = colormap(0.9.*[0,0,1;0,1,0;1,0,0]);

plot(ylist(3,:),ylist(1,:),'k-');

%Add max/min orbit values
hold on
scatter(ylist(3,:),ylist(1,:),15,stab,'filled')
gscatter(rangeI(1,:),rangeI(2,:),orbStab,'gr','__',[5,5],false)
gscatter(rangeI(1,:),rangeI(3,:),orbStab,'gr','++',[5,5],false)

xlabel("beta")
ylabel("I")
xlim([3.5,6]);

title("Bifurcation Diagram")

hold off

%Plot orbit family in beta-T plane
figure(2)

cMap = colormap(0.9.*[0,1,0;1,0,0]);
plot(orbList(3,:),orbList(1,:),'k')
hold on
%Green for stable, red for unstable
scatter(orbList(3,:),orbList(1,:),15,orbStab,'filled')

hold off

xlabel("beta")
ylabel("T")

title("Period of orbits")

%%Plot examples of periodic orbits
figure(3)
sgtitle("Example orbits")

%Select 10 orbits, evenly spaced
orbIndices = floor(linspace(1,nnz(~isnan(orbStab)),10));

i = 1;
while i <= 10
    
    j = orbIndices(i);
    
    %Solve orbit
    [~,~,orb] = MyIVP(@(t,x) rhs(x,orbList(3,j)),[xStar,orbList(2,j)],[0,orbList(1,j)],floor(N/5),'dp45');
    
    %Plot orbit
    subplot(5,2,i)
    
    %Colour by stability
    if(orbStab(j) == 1)
        orbStyle = 'g';
    else
        orbStyle = 'r';
    end
    
    plot(orb(1,:),orb(2,:),orbStyle);
    
    xlim([0,0.35])
    ylim([0,0.75])
    
    xlabel("I")
    ylabel("R")
    
    title("beta = " + num2str(orbList(3,j)));
    
    i = i + 1;
    
end

%%Locate fold point

%Approximately locate fold

foldApprox = zeros(3,1);

i = 2;
while i <= length(orbStab)
    
    if(~isnan(orbStab(i)))
        
        if((orbStab(i) == 1 && orbStab(i-1) == 2) || (orbStab(i) == 2 && orbStab(i-1) == 1))
            foldApprox = orbList(:,i);
            break
        end
        
    end
    
    i = i + 1;
    
end

%For reference
%M = @(t1,x0,p) MyIVP(@(t,x) rhs(x,p),x0,[0,t1],N,'dp45');
%f = @(x) M(x(1),[xStar;x(2)],x(3)) - [xStar;x(2)];
%JM = @(x) MyJacobian(@(y) M(x(1),y,x(3)), [xStar;x(2)], jacobianAccuracy);

%Equations for fold point
res = @(x) [f(x(1:3)); JM(x(1:3)) * x(4:5); x(4:5)' * x(4:5) - 1];
dres = @(x) MyJacobian(res,x,jacobianAccuracy*1e-6);

%Find initial guess for v
[eVec,~] = eigs(JM(foldApprox));

%Initial guess
x0 = [foldApprox;eVec(:,2)];

%Converge to fold point
solFold = MySolve(res,x0,dres,'tol',5e-3);

%Print result
disp(["Fold at beta = ", num2str(solFold(3))])

figure(2)
hold on
plot(solFold(3),solFold(1),'o','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[1,0,1],'MarkerSize',4)

%Legend
leg = zeros(2, 1);
leg(1) = plot(NaN,NaN,'or','MarkerFaceColor','r');
leg(2) = plot(NaN,NaN,'og','MarkerFaceColor','g');
leg(3) = plot(NaN,NaN,'o','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[1,0,1]);
legend(leg, 'Unstable','Stable','Fold','Location','northeast');

hold off

%Find min and max of fold orbit

%Solve for one full orbit
[~,~,orbFold] = MyIVP(@(t,x) rhs(x,solFold(3)),[xStar,solFold(2)],[0,solFold(1)],N,'dp45');

%Calculate min and max
rangeIFold = [0;0];
rangeIFold(1) = min(orbFold(1,:));
rangeIFold(2) = max(orbFold(1,:));

%Plot on bifurcation diagram
figure(1)
hold on

plot(solFold(3),rangeIFold(1),'_','Color',[1,0,1])
plot(solFold(3),rangeIFold(2),'+','Color',[1,0,1])

%Legend
leg = zeros(5, 1);
leg(1) = plot(NaN,NaN,'or','MarkerFaceColor','r');
leg(2) = plot(NaN,NaN,'ob','MarkerFaceColor','b');
leg(3) = plot(NaN,NaN,'og','MarkerFaceColor','g');
leg(4) = plot(NaN,NaN,'k+');
leg(5) = plot(NaN,NaN,'k_');
legend(leg, 'Unstable','Saddle','Stable','Orbit Max','Orbit Min','Location','northwest');

hold off

%%Track fold point - INCOMPLETE
%
% M2 = @(t1,x0,p) MyIVP(@(t,x) rhs2P(x,p(1),(2)),x0,[0,t1],N,'dp45');
%
% g = @(x) M2(x(1),[xStar;x(2)],x(3:4)) - [xStar;x(2)];
% JM2 = @(x) MyJacobian(@(y) M2(x(1),y,x(3:4)), [xStar;x(2)], jacobianAccuracy);
%
% res2 = @(x) [g(x(1:4)); JM2(x(1:4)) * x(5:6); x(5:6)' * x(5:6) - 1];
% dres2 = @(x) MyJacobian(res,x,jacobianAccuracy);
%
% tan1 = [0;0;0;1;0;0];
%
% [eVec,~] = eigs(JM2([foldApprox;gamma0]));
% fInit0 = [foldApprox;gamma0;eVec(:,2)];
% trackFoldList(:,:,1) = MyTrackCurve(res2,fInit0,tan1,'sMax',0.01,'nMax',10,'correctGuess',false);
% trackFoldList(:,:,2) = MyTrackCurve(res2,fInit0,-tan1,'sMax',0.01,'nMax',10,'correctGuess',false);