function [xend, t, xt] = MyIVP(f,x0,tspan,N,METHOD)

if(~(exist("METHOD","var")))
    METHOD = "rk4";
end

if(strcmp(METHOD,"dp45"))
    
    t = linspace(tspan(1),tspan(2),N+1);
    maxIter = length(t);
    h = t(2) - t(1);
    
    %Intialize spatial arrays
    xt = zeros(length(x0),maxIter);
    %Intial Values
    xt(:,1) = x0;
    
    i = 1;
    while (i < maxIter)
        
        k1 = f(t(i),xt(:,i));
        k2 = f(t(i)+(1/5).*h, xt(:,i) + h.*((1/5).*k1));
        k3 = f(t(i)+(3/10).*h, xt(:,i) + h.*((3/40).*k1 + (9/40).*k2));
        k4 = f(t(i)+(4/5).*h, xt(:,i) + h.*((44/45).*k1 + (-56/15).*k2 + (32/9).*k3));
        k5 = f(t(i)+(8/9).*h, xt(:,i) + h.*((19372/6561).*k1 + (-25360/2187).*k2 + (64448/6561).*k3 + (-212/729) .* k4));
        k6 = f(t(i)+h, xt(:,i) + h.*((9017/3168).*k1 + (-355/33).*k2 + (46732/5247).*k3 + (49/176).*k4 + (-5103/18656).*k5));
        
        xt(:,i+1) = xt(:,i) + h.*((35/384).*k1 + 0.*k2 + (500/1113).*k3 + (125/192).*k4 +(-2187/6784).*k5 + (11/84).*k6);
              
        i = i+1;
        
    end

    xend = xt(:,end);
    
else
    
    t = linspace(tspan(1),tspan(2),N+1);
    maxIter = length(t);
    h = t(2) - t(1);
    
    %Intialize spatial arrays
    xt = zeros(length(x0),maxIter);
    %Intial Values
    xt(:,1) = x0;
    
    i = 1;
    while (i<maxIter)
        %k1
        k1 = f(t(i),xt(:,i));
        %k2
        k2 = f(t(i)+(h./2),xt(:,i)+h*(k1./2));
        %k3
        k3 = f(t(i)+(h./2),xt(:,i)+h*(k2./2));
        %k4
        k4 = f(t(i)+h,xt(:,i)+h*k3);
        
        xt(:,i+1) = xt(:,i) + (h/6) .* (k1 + 2*k2 + 2*k3 + k4);
        
        i = i + 1;
        
    end
    
    xend = xt(:,end);
    
end

end