clear

sir_model_680029911;

cMap = colormap(0.9.*[0,0,1;0,1,0;1,0,0]);

jacobianAccuracy = 1e-6;

%Locate the initial equilibrium when B = 6, initial guess I=[0.25;0.25]

B = 6;
f = @(I) rhs(I,B);
df = @(I) MyJacobian(f,I,jacobianAccuracy);

I0 = Solve(f,[0.25;0.25],df);

%Locate equilibria using pseudo arclength continuation

y0 = [I0;B];
yTan0 = [0;0;-1];

g = @(y) rhs(y(1:2),y(3));
dg = @(y)MyJacobian(g,y,jacobianAccuracy);

ylist = MyTrackCurve(g,dg,y0,yTan0,'stop',@(y) (y(3)<3.5 || y(1)<0),'sMax',0.01,'nMax',270);

%%Compute stability; 0 = saddle, 1 = stable, 2 = unstable
stab = NaN(1,length(ylist));

for i = 1:length(ylist)
    
    if(~isnan(ylist(:,i)))
        
        J = dg(ylist(:,i));
        e = eigs(J(1:2,1:2));
        
        if (all(real(e)<0))
            stab(i) = 1;
        elseif (all(real(e)>0))
            stab(i) = 2;
        else
            stab(i) = 0;
        end
        
    end
    
end

plot(ylist(3,:),ylist(1,:),'k-');
xlabel("beta")
ylabel("I")
xlim([3.5,6]);

hold on
scatter(ylist(3,:),ylist(1,:),15,stab,'filled')

leg = zeros(3, 1);
leg(1) = plot(NaN,NaN,'or','MarkerFaceColor','r');
leg(2) = plot(NaN,NaN,'ob','MarkerFaceColor','b');
leg(3) = plot(NaN,NaN,'og','MarkerFaceColor','g');
legend(leg, 'Unstable','Saddle','Stable');

hold off

save('ylist_stab','ylist','stab')