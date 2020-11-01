clear

sir_model_680029911;

%Resolution of jacobian estimate
jacobianAccuracy = 1e-6;

%Plot limits
iRange = [0,0.35];
rRange = [0,0.75];

%Sigma funciton from model
F = @(I) sigm(I,gamma0);

%%%%%%%%%%% Values in the following section are hand selected to provide robust plots and appropriate trajectories/orbits %%%%%%%%%%%

%list of beta values for each plot
%3.7                    Stable
%4.053                  Unstable
%4.14                   Saddle + 2 Unstable
%4.9                    Unstable
%5.5                    Stable, two orbits
%5.8                    Stable, no orbits
betaList = [3.7,4.05,4.14,4.9,5.5,5.8];

%List of initial guesses to locate equilibria for each beta.
guessesList = {...
    [0.012;0.7],...
    [0.01;0],...
    [[0.031;0],[0.018;0],[0.05;0]],...
    [0.083;0.576],...
    [0.091;0],...
    [0.091;0]};

%List of points to plot an orbit from, including period: [I;R;T]
orbitList = {...
    [],...
    [0.08372;0.642;30],...
    [0.08372;0.65;27],...
    [0.08372;0.682;15.15],...
    [[0.08372;0.6723;12.0],[0.08372;0.6877;12.712]],...
    []};

%Extra trajectories (Only used for multiple orbits)

extraList = {...
    [],...
    [],...
    [],...
    [],...
    [0.2;0.2],...
    []};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    %Functions defining model for given beta value
    f = @(I) rhs(I(1:2),beta);
    df = @(I) MyJacobian(f,I,jacobianAccuracy);
    
    %system nullclines
    nulcA = @(I) (1 - I) - ((mu0 + sigma0)./(beta)) - ((eta0.*F(I))./(beta.*I));
    nulcB = @(I) (mu0.*I + eta0.*F(I))./(nu0 + sigma0);
    
    %plot nullclines
    plot(iRange(1):0.001:iRange(2),nulcA(iRange(1):0.001:iRange(2)),'b')
    plot(iRange(1):0.001:iRange(2),nulcB(iRange(1):0.001:iRange(2)),'b')
    
    %List of initial guesses, to be converged to true equilibria
    guesses = cell2mat(guessesList(i));
    
    j = 1;
    while j <= length(guesses(1,:))
        
        %Converge guess to true equilibrium point
        eqPoint = Solve(f,guesses(:,j),df);
        
        [eVecs,eVals] = eigs(df(eqPoint));
        
        %Stability Check
        if (all(real(diag(eVals))>0))
            
            %Unstable
            
            %Plot point
            plot(eqPoint(1),eqPoint(2),'ko')
            
            %Calculate and plot a sample trajectory for this point
            [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+0.01,[0,100],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            drawArrow(t,sol0(1,:),sol0(2,:))
            
        elseif (all(real(diag(eVals))<0))
            
            %Stable
            
            %Plot point
            plot(eqPoint(1),eqPoint(2),'k*')
            
            %Calculate and plot a sample trajectory for this point
            [~,t,sol0]=MyIVP(@(t,x) f(x),eqPoint+[0.1;0],[0,100],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            drawArrow(t,sol0(1,:),sol0(2,:))
            
        else
            
            %Saddle
            
            %Plot point
            plot(eqPoint(1),eqPoint(2),'kx')
            
            %Calculate and plot separatrices
            [~,~,separatrix1]=MyIVP(@(t,x) f(x),eqPoint+eVecs(:,1).*1e-6,[0,100],2500,'dp45'); %Unstable
            [~,~,separatrix2]=MyIVP(@(t,x) f(x),eqPoint-eVecs(:,1).*1e-6,[0,100],2500,'dp45');
            [~,~,separatrix3]=MyIVP(@(t,x) f(x),eqPoint+eVecs(:,2).*1e-2,[100,0],2500,'dp45'); %Stable
            [~,~,separatrix4]=MyIVP(@(t,x) f(x),eqPoint-eVecs(:,2).*1e-2,[100,0],2500,'dp45');
            
            plot(separatrix1(1,:),separatrix1(2,:),'r');
            plot(separatrix2(1,:),separatrix2(2,:),'r');
            plot(separatrix3(1,:),separatrix3(2,:),'g');
            plot(separatrix4(1,:),separatrix4(2,:),'g');
            
        end
        
        j = j + 1;
        
    end
    
    %List of start points for additional trajectories;
    extraStarts = cell2mat(extraList(i));
    
    j = 1;
    if(~isempty(extraStarts))
        while j <= length(extraStarts(1,:))
            
            %Simulate trajectory and plot
            [~,t,sol0]=MyIVP(@(t,x) f(x),extraStarts(1:2,j),[0,100],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'k')
            drawArrow(t,sol0(1,:),sol0(2,:))
            
            j = j + 1;
            
        end
    end
    
    %List of points approximately located on an orbit
    orbitGuesses = cell2mat(orbitList(i));
    
    j = 1;
    if(~isempty(orbitGuesses))
        while j <= length(orbitGuesses(1,:))
            
            %Simulate for one orbit and plot
            [~,t,sol0]=MyIVP(@(t,x) f(x),orbitGuesses(1:2,j),[0,orbitGuesses(3,j)],2500,'dp45');
            plot(sol0(1,:),sol0(2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
            
            j = j + 1;
            
        end
    end
    
    %Add legend
    leg = zeros(8, 1);
    leg(1) = plot(NaN,NaN,'k*');
    leg(2) = plot(NaN,NaN,'ko');
    leg(3) = plot(NaN,NaN,'kx');
    leg(4) = plot(NaN,NaN,'k_');
    leg(5) = plot(NaN,NaN,'_','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    leg(6) = plot(NaN,NaN,'b_');
    leg(7) = plot(NaN,NaN,'r_');
    leg(8) = plot(NaN,NaN,'g_');
    legend(leg, 'Sink','Source','Saddle','Trajectory','Orbit','Nullcline','Seperatrix (un)','Seperatrix (stab)','Location','southeast');
    
    
    hold off
    
    i = i + 1;
    
end