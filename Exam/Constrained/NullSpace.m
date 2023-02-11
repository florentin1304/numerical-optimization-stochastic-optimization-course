function [x_star, fxstar, lambda_star, v_star, KKT_gradL_norm, KKT_eq_norm] = ...
    NullSpace(Q, c, A, b, x2)


% INPUTS: 
% Q = symmetric positive semidefinite matrix n by n
% c = n-dimensional vector of the quadratic loss function
% A = equality constraints matrix, size K x n
% b = equality constraints vector, K-dimensional 
% x2 = column vector of size (n-K)


% OUTPUTS:
% xstar = solution 
% fxstar = value of the loss function in xstar
% lambda_star = lagrangian multiplier computed by the function
% v_star = solution of the substitution variable
% KKT_gradL_norm, KKT_eqnorm: norm of the respective corresponding KKT conditions


[m, n] = size(A);
% x_hat special solution of Ax=b --> divide A in A1, A2
A1 = A(:, 1:m);
A2 = A(:, m+1:end);
% find Z such that A*Z=0
Z = [-A1\A2; speye(n-m)];


% xhat initialization
x_hat = [A1\(b - A2 * x2); x2];

% compute v_star as solution of:
% (Z' Q Z)v = -Z'(Q xhat + c )

lhs = Z' * Q * Z;
rhs = -Z' * (Q * x_hat + c);
v_star = lhs\rhs;

% compute xstar
x_star = Z * v_star + x_hat;

% compute lambda_star
lambda_star = (A * A')\(-A *(c + Q * x_star));

% compute fxstar
fxstar = 0.5 * x_star' * Q * x_star + c' * x_star;

KKT_gradL_norm = norm(Q * x_star + c + A' * lambda_star);
KKT_eq_norm = norm(A * x_star - b);

end