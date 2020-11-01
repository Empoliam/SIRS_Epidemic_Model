clear

sir_model_680029911;

jacobianAccuracy = 1e-6;
iRange = [0,0.35];
rRange = [0,0.75];

%F = @(I) gamma0 .* (((2)./(1+exp(-2.*I./gamma0)))-1);
F = @(I) sigm(I,gamma0);

%list of beta values for each plot
%3.7                    Stable
%4.053                  Unstable
%4.14                   Saddle + 2 Unstable
%4.9                    Unstable
%5.5                    Stable

betaList = [3.7,4.05,4.14,4.9,5.5,5.8];

%List of initial guesses to locate equilibria for each beta.
guessesList = {...
    [0.012;0.7],...
    [0.01;0],...
    [[0.031;0],[0.018;0],[0.05;0]],...
    [0.083;0.576],...
    [0.091;0],...
    [0.091;0]};

%loop counter
i = 1;
for beta = betaList
         
    figure(i)
    
    xlabel('I')
    ylabel('R')
    
    xlim(iRange)
    ylim(rRange)
    
    pbaspect([1 1 1])
    
    title(['beta = ', num2str(beta)])
    
    hold on
    
    f = @(I) rhs(I(1:2),beta);
    df = @(I) MyJacobian(f,I,jacobianAccuracy);
        
    nulcA = @(I) (1 - I) - ((mu0 + sigma0)./(beta)) - ((eta0.*F(I))./(beta.*I));
    nulcB = @(I) (mu0.*I + eta0.*F(I))./(nu0 + sigma0);
    
    %plot nullclines
    plot(iRange(1):0.001:iRange(2),nulcA(iRange(1):0.001:iRange(2)),'g')
    plot(iRange(1):0.001:iRange(2),nulcB(iRange(1):0.001:iRange(2)),'g')
    
    guesses = cell2mat(guessesList(i));
    
    j = 1;
    while j <= length(guesses(1,:))
        
        eqPoint = Solve(f,guesses(:,j),df);
        
        [eVecs,eVals] = eigs(df(eqPoint));
        
        %Stability Check
        if (all(real(diag(eVals))>0))
            
            %Unstable
            
            plot(eqPoint(1),eqPoint(2),'ko')
            [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+0.01,[0,100],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            drawArrow(t,sol0(1,:),sol0(2,:))
            
        elseif (all(real(diag(eVals))<0))
            
            %Stable
            
            % interior orbit ~ +[0.102675;0]
            
            plot(eqPoint(1),eqPoint(2),'k*')
            [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+[0.3;0],[0,100],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            drawArrow(t,sol0(1,:),sol0(2,:))
            
%             [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+[0.11;0],[0,100],2500,'dp45');
%             plot(sol0(1,:),sol0(2,:),'c')
%             drawArrow(t,sol0(1,:),sol0(2,:),'lineSpec','c')
%             [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+[0.2;0],[0,100],2500,'dp45');
%             plot(sol0(1,:),sol0(2,:),'b')
%             drawArrow(t,sol0(1,:),sol0(2,:),'lineSpec','b')
%             [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+[0.1;0],[0,100],2500,'dp45');
%             plot(sol0(1,:),sol0(2,:),'r')
%             drawArrow(t,sol0(1,:),sol0(2,:),'lineSpec','r')
            
        else
            
            %Saddle
            
            plot(eqPoint(1),eqPoint(2),'kx')
            
            [~,~,separatrix1]=MyIVP(@(t,x) f(x),eqPoint+eVecs(:,1).*1e-6,[0,100],2500,'dp45');
            [~,~,separatrix2]=MyIVP(@(t,x) f(x),eqPoint-eVecs(:,1).*1e-6,[0,100],2500,'dp45');
            [~,~,separatrix3]=MyIVP(@(t,x) f(x),eqPoint+eVecs(:,2).*1e-2,[100,0],2500,'dp45');
            [~,~,separatrix4]=MyIVP(@(t,x) f(x),eqPoint-eVecs(:,2).*1e-2,[100,0],2500,'dp45');
            
            plot(separatrix1(1,:),separatrix1(2,:),'r');
            plot(separatrix2(1,:),separatrix2(2,:),'r');
            plot(separatrix3(1,:),separatrix3(2,:),'b');
            plot(separatrix4(1,:),separatrix4(2,:),'b');
            
        end
        
        j = j + 1;
        
    end
       
    hold off
    
    i = i + 1;
    
end