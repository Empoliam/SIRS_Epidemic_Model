clear

sir_model_680029911;
cMap = colormap(0.8.*[0,0,1;0,1,0;1,0,0]);
jacobianAccuracy = 1e-6;

load('ylist_stab.mat');

%Equations defining model, accepts one parameter
f = @(I) rhs(I(1:2),I(3));
df = @(I) MyJacobian(@(a) rhs(a,I(3)),I(1:2),jacobianAccuracy);

%find approximate bifurcation locations using MyTrackCurve results
hopfApprox = []; %Array of approximate hopf locations
foldApprox = []; %Array of approximate fold locations
i = 2;
while i <= length(stab)
    
    if(~isnan(stab(i)))
        
        if((stab(i) == 1 && stab(i-1) == 2) || (stab(i) == 2 && stab(i-1) == 1)) %Searching for locations where a stable node becomes unstable
            hopfApprox(end + 1) = i;
        elseif((stab(i) == 0 && stab(i-1) ~= 0) || (stab(i) ~= 0 && stab(i-1) == 0)) %Searching for locations where a saddle node becomes unstable
            foldApprox(end + 1) = i;
        end
        
    end
    
    i = i + 1;
    
end

%System of equations defining a fold point
ffold = @(I) [f(I(1:3)); df(I(1:3)) * I(4:5); I(4:5)' * I(4:5) - 1];
dffold = @(I) MyJacobian(ffold,I,jacobianAccuracy);

foldList = []; %Array of true fold locations

for i = foldApprox
         
    [eVec,~] = eigs(df(ylist(:,i))); %Calculate eigenvectors of jacobian
    I0 = [ylist(:,i);eVec(:,2)]; %Concatenate initial guess with eigenvector to create initial guess for fold point system
    sol = Solve(ffold,I0,dffold);
    foldList = [foldList, sol(1:3)]; %Add to list of solutions
    
end

%System of equations defining a hopf point
fhopf = @(I) [f(I);trace(df(I))];
dfhopf = @(I) MyJacobian(fhopf,I,jacobianAccuracy);

hopfList = []; %Array of true hopf locations

for i = hopfApprox
    
    I0 = ylist(:,i);
    hopfList = [hopfList , Solve(fhopf,I0,dfhopf)]; %Solve and add to list of solutions
    
end

%Plot bifurcation diagram
plot(ylist(3,:),ylist(1,:),'c-');

hold on

scatter(ylist(3,:),ylist(1,:),10,stab,'filled') %Add stabilities
scatter(hopfList(3,:),hopfList(1,:),20,[1,1,0],'filled') %Add hopf points
scatter(foldList(3,:),foldList(1,:),20,[1,0,1],'filled') %Add fold points

%Legend
leg = zeros(5, 1);
leg(1) = plot(NaN,NaN,'or','MarkerFaceColor','r');
leg(2) = plot(NaN,NaN,'ob','MarkerFaceColor','b');
leg(3) = plot(NaN,NaN,'og','MarkerFaceColor','g');
leg(4) = plot(NaN,NaN,'o','MarkerFaceColor',[1,1,0],'MarkerEdgeColor',[1,1,0]);
leg(5) = plot(NaN,NaN,'o','MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[1,0,1]);
legend(leg, 'Unstable','Saddle','Stable','Hopf','Fold');

hold off

xlabel("beta")
ylabel("I")

xlim([3.5,6]);

%%3b

%Hopf bifurcation tracking

stopCond = @(y) (y(4) > 0.1 || y(4) < 0 || y(3) > 12 || y(3) < 0); %Condition for stoppping curve tracking

%Equations defining model, accepting two parameters
g = @(y) rhs2P(y(1:2),y(3),y(4));
dg = @(I) MyJacobian(@(a) rhs2P(a,I(3),I(4)),I(1:2),jacobianAccuracy);

%New equations for hopf points
ghopf = @(I) [g(I);trace(dg(I))];
dghopf = @(I) MyJacobian(ghopf,I,jacobianAccuracy);

%Initial tangent guess, parameters decreasing
tan0 = [0;0;-1;-1];
tan0 = tan0/norm(tan0);

%Track hopf bifurcation from one hopf point
hInit0 = [hopfList(:,1);gamma0];
trackHopfList(:,:,1) = MyTrackCurve(ghopf,hInit0,tan0,'stop',stopCond,'sMax',0.05,'nMax',500,'correctGuess',true);
trackHopfList(:,:,2) = MyTrackCurve(ghopf,hInit0,-tan0,'stop',stopCond,'sMax',0.05,'nMax',500,'correctGuess',true);

%Track hopf bifurcation from other hopf point
hInit1 = [hopfList(:,2);gamma0];
trackHopfList(:,:,3) = MyTrackCurve(ghopf,hInit1,tan0,'stop',stopCond,'sMax',0.05,'nMax',500,'correctGuess',true);
trackHopfList(:,:,4) = MyTrackCurve(ghopf,hInit1,-tan0,'stop',stopCond,'sMax',0.05,'nMax',500,'correctGuess',true);

%Check that hopf points are not saddles
i = 1;
while i < 4
    
    j = 1;
    while j < length(trackHopfList(1,:,i))
        
        pointCheck = trackHopfList(:,j,i);
        
        if(~isnan(pointCheck))
            
            eVals = eigs(dg(pointCheck)); %Find eigenvalues of jacobian at point
            
            if(all(imag(eVals) == 0))
                trackHopfList(:,j,i) = NaN(4,1); %Discard points with only real eigenvalues (saddles)
            end
            
        end
        
        j = j + 1;
    end
    
    i = i + 1;
    
end

%%Fold bifurcation tracking

% Equations defining fold bifurcations, accepting two parameters
gfold = @(I) [g(I(1:4)); dg(I(1:4)) * I(5:6); I(5:6)' * I(5:6) - 1];
dgfold = @(I) MyJacobian(gfold,I,jacobianAccuracy);

%Initial tangent guess, parameters decreasing
tan1 = [0;0;-1;-1;0;0];
tan1 = tan1/norm(tan1);

%Track from first fold point
[vec0,~] = eigs(dg([foldList(:,1);gamma0])); %Find eigenvectors at point, for initial v guess
fInit0 = [foldList(:,1);gamma0;vec0(:,2)]; %Concatenate starting point with v guess
trackFoldList(:,:,1) = MyTrackCurve(gfold,fInit0,tan1,'stop',stopCond,'sMax',0.01,'nMax',350,'correctGuess',true);
trackFoldList(:,:,2) = MyTrackCurve(gfold,fInit0,-tan1,'stop',stopCond,'sMax',0.01,'nMax',350,'correctGuess',false);

%Track from other fold point
[vec1,~] = eigs(dg([foldList(:,1);gamma0])); %Find eigenvectors at point, for initial v guess
fInit1 = [foldList(:,2);gamma0;vec1(:,2)]; %Concatenate starting point with v guess
trackFoldList(:,:,3) = MyTrackCurve(gfold,fInit1,tan1,'stop',stopCond,'sMax',0.01,'nMax',350,'correctGuess',true);
trackFoldList(:,:,4) = MyTrackCurve(gfold,fInit1,-tan1,'stop',stopCond,'sMax',0.01,'nMax',350,'correctGuess',true);

%Plot bifurcation curves
figure()

plot(trackHopfList(3,:),trackHopfList(4,:),'r')
hold on
plot(trackFoldList(3,:),trackFoldList(4,:),'b')

%Legend
leg = zeros(2, 1);
leg(1) = plot(NaN,NaN,'-r');
leg(2) = plot(NaN,NaN,'-b');
legend(leg, 'Hopf','Fold');
hold off

xlabel('beta')
xlim([0,12])
ylabel('gamma')
ylim([0,0.07])


save('hopf_fold','hopfList','foldList')
