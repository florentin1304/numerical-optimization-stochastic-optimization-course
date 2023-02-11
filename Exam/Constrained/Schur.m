function [xstar, fxstar, lambda_star, KKT_gradL_norm, KKT_eq_norm] = ...
    Schur(Q, c, A, b)


% INPUTS: 
% Q = symmetric positive semidefinite matrix n by n
% c = n-dimensional vector of the quadratic loss function
% A = equality constraints matrix, size K x n
% b = equality constraints vector, K-dimensional 


% OUTPUTS:
% xstar = solution 
% fxstar = value of the loss function in xstar
% lambda_star = lagrangian multiplier computed by the function
% KKT_gradL_norm, KKT_eqnorm: norm of the respective corresponding KKT conditions

[n, ~] = size(Q);
Q_inv = Q\eye(n);

% Schur complement --> Q_hat is symmetric positive def under assumptions
Q_hat = A * Q_inv * A';


% find lambda_star solving the linear system:
% Schur * lambda = -b -A * Qinv * c

rhs = -b - A * Q_inv * c;  % right-hand-side of the prev equation
lambda_star = Q_hat\rhs;

% xstar given lambda_star
xstar = Q_inv * (-c - A' * lambda_star);


% fxstar
fxstar = 0.5 * xstar' * Q * xstar + c' * xstar;

KKT_gradL_norm = norm(Q * xstar + c + A' * lambda_star);
KKT_eq_norm = norm(A * xstar - b);

end