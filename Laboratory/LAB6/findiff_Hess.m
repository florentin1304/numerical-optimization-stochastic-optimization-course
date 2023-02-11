function [hess_fx] = findiff_Hess(f,x,h)
%FINDIFF_HESS Approximation of f's hessian in x

    n = length(x);
    E = eye(n,n);
    hess_fx = zeros(n,n);
    fx = f(x);
    for i=1:n
        % from i to n ?
        for j=1:n
            if i == j
                hess_fx(i,j) =  (f(x+h*E(:,i))  - (2*fx)+ f(x-h*E(:,i)) ) / (h^2);
            else
                hess_fx(i,j) = ( f( x + h*E(:,i) + h*E(:,j) ) - f(x+h*E(:,i)) - f(x+h*E(:,j)) + fx ) / (h^2);
            end
        end
    end
    
end

