clear

cMap = colormap(0.9.*[0,0,1;0,1,0;1,0,0]);

sir_model_680029911;

F = @(I) gamma0 .* (((2)./(1+exp(-2.*I./gamma0)))-1);

%list of beta values for each plot
%4.14 Saddle + 2 Unstable
%5.5 Stable
%4.7 Unstable
betaList = [4.14,5.5,4.7];

%List of initial guesses to locate equilibria for each beta.
guessesList = {...
    [[0.031;0],[0.018;0],[0.05;0]],...
    [0.091;0],...
    [0.091;0]};

%loop counter
i = 1;
for beta = betaList
        
    figure(i)
    
    hold on
    
    f = @(I) rhs(I(1:2),beta);
    df = @(I) MyJacobian(f,I,1e-6);
    
    nulcA = @(I) (1 - I) - ((mu0 + sigma0)./(beta)) - ((eta0.*F(I))./(beta.*I));
    nulcB = @(I) (mu0.*I + eta0.*F(I))./(nu0 + sigma0);
    
    %plot nullclines
    plot(0:0.001:0.3,nulcA(0:0.001:0.3),'g')
    plot(0:0.001:0.3,nulcB(0:0.001:0.3),'g')
    
    guesses = cell2mat(guessesList(i));
    
    j = 1;
    while j <= length(guesses(1,:))
        
        eqPoint = Solve(f,guesses(:,j),df);
        
        [eVecs,eVals] = eigs(df(eqPoint));
        
        %Stability Check
        if (all(diag(eVals)>0))
            
            %Unstable
            
            plot(eqPoint(1),eqPoint(2),'o')
            [~,~,sol0]=MyIVP(@(t,x) f(x),eqPoint+0.01,[0,50],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            
        elseif (all(diag(eVals)<0))
            
            %Stable
            
            plot(eqPoint(1),eqPoint(2),'*')
            [~,~,sol0]=MyIVP(@(t,x) f(x),eqPoint+0.01,[0,50],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            
            
        else
            
            %Saddle
            
            plot(eqPoint(1),eqPoint(2),'x')
            
            [~,~,separatrix1]=MyIVP(@(t,x) f(x),eqPoint+eVecs(:,1).*1e-6,[0,50],2500,'dp45');
            [~,~,separatrix2]=MyIVP(@(t,x) f(x),eqPoint-eVecs(:,1).*1e-6,[0,50],2500,'dp45');
            [~,~,separatrix3]=MyIVP(@(t,x) f(x),eqPoint+eVecs(:,2).*1e-2,[100,0],2500,'dp45');
            [~,~,separatrix4]=MyIVP(@(t,x) f(x),eqPoint-eVecs(:,2).*1e-2,[100,0],2500,'dp45');
            
            plot(separatrix1(1,:),separatrix1(2,:),'r');
            plot(separatrix2(1,:),separatrix2(2,:),'r');
            plot(separatrix3(1,:),separatrix3(2,:),'b');
            plot(separatrix4(1,:),separatrix4(2,:),'b');
            
        end
        
        j = j + 1;
        
    end
    
    xlabel('I')
    ylabel('R')
    
    xlim([0,0.3])
    ylim([0,0.7])
    
    hold off
    
    i = i + 1;
    
end