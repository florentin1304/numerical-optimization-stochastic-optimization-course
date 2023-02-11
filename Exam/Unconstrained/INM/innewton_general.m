function [xk, fk, gradfk_norm, k, xseq, btseq, fseq, gfseq] = ...
    innewton_general(x0, f, gradf, Hessf, kmax, ...
    tolgrad, toldiffgrad, c1, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit)

% INPUTS:
% x0 = column vector of dimension nx1;
% f = function handle;
% gradf = function handle that describes the gradient of f;
% Hessf = function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = stopping criterion;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy;
% FDgrad = 'fw', 'c' 
% FDHess = 'fw', 'Jfw', 'Jc', 'MF'
% h = approximation step for FD;
% fterms = f. handle "@(gradfk, k) ..." that returns the forcing term
% eta_k at each iteration
% pcg_maxit = maximum number of iterations for the pcg solver.


% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.


switch FDgrad
    case 'fw'
        gradf = @(x) findiff_grad(f, x, h, 'fw');
        
    case 'c'
        gradf = @(x) findiff_grad(f, x, h, 'c');
        
    otherwise
        % WE USE THE INPUT FUNCTION HANDLE gradf...
        %
        % THEN WE DO NOT NEED TO WRITE ANYTHING!
end

if isequal(FDgrad, 'fw') || isequal(FDgrad, 'c')
    switch FDHess
        case 'fw'
            Hessf = @(x) findiff_Hess(f, x, sqrt(h));
        case 'MF'
            Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
        otherwise
    end
else
    switch FDHess
        case 'fw'
            Hessf = @(x) findiff_Hess(f, x, sqrt(h));
        case 'Jfw'
            Hessf = @(x) findiff_J(gradf, x, h, 'fw');

        case 'Jc'
            Hessf = @(x) findiff_J(gradf, x, h, 'c');

        case 'MF'
            Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
        otherwise
    end
end



% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk)  fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
fseq = zeros(1, kmax);
gfseq = zeros(1, kmax);

same = 0;
n = size(x0);
xk = x0;
fk = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);

fseq(1) = fk;
gfseq(1) = gradfk_norm;

while k < kmax && gradfk_norm >= tolgrad
    % Hessf(xk) p = - graf(xk)
    
    eta_k = fterms(gradfk, k);

    switch FDHess
        case 'MF'
            Hessfk_pk = @(p) Hessf_pk(xk, p);

            pk = pcg(Hessfk_pk, -gradfk, eta_k, pcg_maxit);

        otherwise
            Hessfk = Hessf(xk);
            if(rank(Hessfk) < n)
                error("Hessian not singular !")
            end
            e = eig(Hessfk);
            
            if all(e > 0)
                disp("Hessian PD")
            else
                disp("Hessian non-PD")
            end

            c = cond(Hessfk);
            fprintf("Cond(H): %f\n", c);

            pk = pcg(Hessf(xk), -gradfk, eta_k, pcg_maxit);
    end
    
    
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, xk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm_old = gradfk_norm;
    gradfk_norm = norm(gradfk);
    fprintf("Iteration %d gf_norm: %f\n", k, gradfk_norm);

    % Increase the step by one
    k = k + 1;
    
    xseq(:, k) = xk;
    btseq(k) = bt;
    fseq(k) = fk;
    gfseq(k) = gradfk_norm;


    if abs(gradfk_norm_old - gradfk_norm) < toldiffgrad
        same = same + 1;
        if same > 4
            break;          
        end
    else
        same = 0;
    end
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
fseq = fseq(1:k);
gfseq = gfseq(1:k);
end