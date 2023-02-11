function[xk, fk, gradfk_norm, k, xseq, btseq, fseq, gfseq] = steepest_desc_bt(x0, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax)

% INPUTS:
% x0 = column vector of size nx1;
% f = function handle that describes the function;
% gradf =  function handle that describes the gradient of the function;
% alpha0 =  initial factor;
% kmax = maximum number of iterations; 
% tolgrad = stopping criterion for the norm of the gradient;
% c1 = factor for Armijo condition;
% rho = factor used for reducing alpha;
% btmax = max number of steps related to the update of alpha in the
% backtrack;


% OUTPUTS:
% xk = last value computed by the function;
% fk = value of f in xk; 
% gradfk_norm = value of the norm of gradf in xk;
% k = index of last iteration performed;
% xseq = matrix of size nxk which columns are xk computed;
% btseq = vector of size 1xk where elements are the number of backtracking
% iterations at each step of optimization;

xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
fseq = zeros(1, kmax);
gfseq = zeros(1, kmax);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradf(xk));

fseq(1) = fk;
gfseq(1) = gradfk_norm;

farmijo = @(fk, alpha, gradfk, pk)...
            fk + c1 * alpha * (gradfk)' * pk;

while k < kmax && gradfk_norm >= tolgrad
    % descent direction
    pk = -gradf(xk);
    alpha = alpha0;
    xnew = xk + alpha * pk;
    fnew = f(xnew);

    % backtrack
    bt = 0;
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk)
        % reduce alpha
        alpha = rho * alpha;
        xnew = xk + alpha * pk;
        fnew = f(xnew);

        bt = bt + 1;
    end

    % uodate with new values
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    k = k + 1;
    
    xseq(:, k) = xk;
    btseq(k) = bt; 
    fseq(k) = fk;
    gfseq(k) = gradfk_norm;

    fprintf("Iteration %d gf_norm: %f\n", k, gradfk_norm);
end

xseq = xseq(:, 1:k);
btseq = btseq(1:k);
fseq = fseq(1:k);
gfseq = gfseq(1:k);
end