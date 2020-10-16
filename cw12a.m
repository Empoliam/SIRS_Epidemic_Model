clear

sir_model_680029911;

F = @(I) gamma0 .* (((2)./(1+exp(-2.*I./gamma0)))-1);

%%saddle region

beta0 = 4.14;
f = @(I) rhs(I(1:2),beta0);
df = @(I) MyJacobian(f,I,1e-6);

nulcA = @(I) (1 - I) - ((mu0 + sigma0)./(beta0)) - ((eta0.*F(I))./(beta0.*I));
nulcB = @(I) (mu0.*I + eta0.*F(I))./(nu0 + sigma0);

E0(:,1) = Solve(f,[0.031;0],df);
E0(:,2) = Solve(f,[0.018;0],df);
E0(:,3) = Solve(f,[0.05;0],df);

[saddleV,saddleE] = eigs(df(E0(:,1)));
eigs(df(E0(:,2)));
eigs(df(E0(:,3)));

[~,~,sol0]=MyIVP(@(t,x) f(x),E0(:,1)+0.001,[0,50],2500,'dp45');
[~,~,sol1]=MyIVP(@(t,x) f(x),E0(:,2)+0.001,[0,50],2500,'dp45');
[~,~,sol2]=MyIVP(@(t,x) f(x),E0(:,3)+0.001,[0,50],2500,'dp45');

[~,~,separatrix1]=MyIVP(@(t,x) f(x),E0(:,1)+saddleV(:,1).*1e-6,[0,50],2500,'dp45');
[~,~,separatrix2]=MyIVP(@(t,x) f(x),E0(:,1)-saddleV(:,1).*1e-6,[0,50],2500,'dp45');
[~,~,separatrix3]=MyIVP(@(t,x) f(x),E0(:,1)+saddleV(:,2).*1e-2,[100,0],25000,'dp45');
[~,~,separatrix4]=MyIVP(@(t,x) f(x),E0(:,1)-saddleV(:,2).*1e-2,[100,0],25000,'dp45');

plot(E0(1,:),E0(2,:),'x')
hold on
%plot(sol0(1,:),sol0(2,:),'k')
%plot(sol1(1,:),sol1(2,:),'k')
%plot(sol2(1,:),sol2(2,:),'k')

plot(0:0.001:0.3,nulcA(0:0.001:0.3),'g')
plot(0:0.001:0.3,nulcB(0:0.001:0.3),'g')

plot(separatrix1(1,:),separatrix1(2,:),'r');
plot(separatrix2(1,:),separatrix2(2,:),'r');
plot(separatrix3(1,:),separatrix3(2,:),'b');
plot(separatrix4(1,:),separatrix4(2,:),'b');

xlabel('I')
ylabel('R')

xlim([0,0.3])
ylim([0,0.7])

hold off


%%stable region
figure()

beta0 = 5.5;
f = @(I) rhs(I(1:2),beta0);
df = @(I) MyJacobian(f,I,1e-6);

nulcA = @(I) (1 - I) - ((mu0 + sigma0)./(beta0)) - ((eta0.*F(I))./(beta0.*I));
nulcB = @(I) (mu0.*I + eta0.*F(I))./(nu0 + sigma0);

E0 = zeros(0,0);
E0(:,1) = Solve(f,[0.091;0],df);
eigs(df(E0(:,1)));

init = E0(:,1)+0.01;

[~,~,sol0]=MyIVP(@(t,x) f(x),init,[0,50],2500,'dp45');
plot(E0(1,:),E0(2,:),'x')
hold on
plot(sol0(1,:),sol0(2,:),'k')

plot(0:0.001:0.3,nulcA(0:0.001:0.3),'g')
plot(0:0.001:0.3,nulcB(0:0.001:0.3),'g')

xlabel('I')
ylabel('R')

hold off

%%unstable region
figure()

beta0 = 4.7;
f = @(I) rhs(I(1:2),beta0);
df = @(I) MyJacobian(f,I,1e-6);

nulcA = @(I) (1 - I) - ((mu0 + sigma0)./(beta0)) - ((eta0.*F(I))./(beta0.*I));
nulcB = @(I) (mu0.*I + eta0.*F(I))./(nu0 + sigma0);

E0 = zeros(0,0);
E0(:,1) = Solve(f,[0.091;0],df);
eigs(df(E0(:,1)));

init = E0(:,1)+0.01;

[~,~,sol0]=MyIVP(@(t,x) f(x),init,[0,50],2500,'dp45');
plot(E0(1,:),E0(2,:),'x')
hold on
plot(sol0(1,:),sol0(2,:),'k')

plot(0:0.001:0.3,nulcA(0:0.001:0.3),'g')
plot(0:0.001:0.3,nulcB(0:0.001:0.3),'g')

xlabel('I')
ylabel('R')

hold off