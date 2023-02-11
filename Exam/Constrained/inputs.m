function [Q, c, A, b] = inputs(n, K)

v = -1*ones(1, n);
v2 = 2*ones(1, n);
c = ones(n, 1);
b = ones(K, 1);


A = zeros(K, n);

for i=1:K
    v1 = zeros(n, 1);
    v1(i:K:end) = 1;
    A(i, :) = v1';
end    

if n>1e4
    e = ones(n,1);
    Q = spdiags([-e 2*e -e],-1:1,n,n);
else
    dg = diag(v, -1) + diag(v, 1);
    Q = diag(v2) + dg(1:n, 1:n);
end


