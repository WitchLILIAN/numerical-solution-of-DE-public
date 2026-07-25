function [x,U] = two_point_reaction(a,b,ua,ub,f,n,K)

if nargin < 7
    K = 0;
end

h = (b-a)/n;
h2 = h*h;

A = sparse(n-1,n-1);
F = zeros(n-1,1);
x = zeros(n-1,1);

for i = 1:n-2
    A(i,i) = -2/h2 - K;
    A(i+1,i) = 1/h2;
    A(i,i+1) = 1/h2;
end
A(n-1,n-1) = -2/h2 - K;

for i = 1:n-1
    x(i) = a + i*h;
    F(i) = feval(f,x(i));
end

F(1) = F(1) - ua/h2;
F(n-1) = F(n-1) - ub/h2;

U = A\F;

return
