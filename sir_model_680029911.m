%% 
% Model script for CW1 of MTH3039 autumn 2020
% vars I,R (infected, recovered)
% free pars beta
% f=rhs(x,beta) where x(1)=I, x(2)=R
clear
format compact
%% Set your personal parameter value for EL
% insert your student number xxxyyyzzz (usually starts with a 6)
% rename as spiking_model_xxxyyyzzz.m and call this in your scripts for
% each questions
SNumber=680029911;

if SNumber==0
   error(['Enter your 9 digit student number in your local copy of spiking_model.m   ',...
            'See instruction sheet or contact James Rankin at j.a.rankin@exeter.ac.uk if unsure.']); 
end
rng(SNumber)
sigma0=0.04+rand(1)*0.01;
disp('Your personal value of parameter sigma:');
format longg
disp(sigma0)
format short

%% Fixed parameters (suffixed with 0)
beta0=6;
mu0=1;
nu0=0.2;
eta0=2.5;
gamma0=0.0225;

S=@(I,R)1-I-R;% equation for S
sigm=@(I,gamma)gamma.*(2./(1+exp(-2./gamma.*I))-1); % smooth saturation function
sirs=@(I,R,beta,mu,sigma,nu,eta,gamma)[beta*I.*S(I,R)-(mu+sigma)*I-eta*sigm(I,gamma);...
    mu*I+eta*sigm(I,gamma)-(nu+sigma)*R];

% put the right hand side in a standard form R^2[x] x R[p] -> R^2[f(x,p)]
rhs=@(x,p)sirs(x(1,:),x(2,:),p,mu0,sigma0,nu0,eta0,gamma0);

% put the right hand side in a standard form R^2[x] x R^2[p1,p2] -> R^2[f(x,p1,p2)]
rhs2P=@(x,p1,p2)sirs(x(1,:),x(2,:),p1,mu0,sigma0,nu0,eta0,p2);
    


