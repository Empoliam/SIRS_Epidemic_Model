function o = MyJacobian(f,x,h)

kMax = length(x);
nf = length(f(zeros(kMax,1)));
I = eye(kMax);

o = zeros(nf,kMax);

for k = 1:kMax
    o(:,k) = (f(x + I(:,k).*h) - f(x - I(:,k).*h))/(2.*h);
end

end