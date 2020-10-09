clear

sir_model_680029911;

%Locate the initial equilibrium when B = 6, initial guess I=[0.25;0.25]

B = 6;
f = @(I) rhs(I,B);
df = @(I) MyJacobian(f,I,1e-6);

I0 = Solve(f,[0.25;0.25],df);

%Locate equilibria using pseudo arclength continuation

y0 = [I0;B];

g = @(y) rhs(y(1:2),y(3));
dg = @(y)MyJacobian(f,y,1e-6);

ylist = MyTrackCurve(g,dg,y0,[0;0;-1],'stop',@(y) (y(3)<3.5 || y(1)<0),'sMax',0.01,'nMax',270);

plot(ylist(3,:),ylist(1,:),'k-');
xlabel("beta")
ylabel("I")
xlim([3.5,6]);


%%Compute stability; 0 = saddle, 1 = stable, 2 = unstable
stab = zeros(1,length(ylist));

for i = 1:length(ylist)
    
    e = real(eigs(df(ylist(1:2,i))));
    
    if (e < 0) 
        stab(i) = 1;
    elseif (e > 0)
        stab(i) = 2;
    else
        stab(i) = 0;
    end
    
end

cMap = colormap(0.9.*[1,0.5,0;0,1,0;1,0,0]);

hold on
scatter(ylist(3,:),ylist(1,:),10,stab,'filled')
hold off