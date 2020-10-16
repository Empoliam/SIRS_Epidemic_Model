clear

sir_model_680029911;

f = @(I) rhs(I(1:2),5.6);
df = @(I) MyJacobian(f,I,1e-6);

A = [1,0;0,1];
B = [1,0;0,-1];
C = [-1,0;0,-1];

M = B;


if (all(diag(M)>0))
    2
elseif (all(diag(M)<0))
    1
else
    0
end