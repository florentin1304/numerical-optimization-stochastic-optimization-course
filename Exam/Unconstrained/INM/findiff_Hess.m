function [Hessfx] = findiff_Hess(f, x, h)

% INPUTS:
% f = function handle that describes a function R^n->R;
% x = n-dimensional column vector;
% h = the h used for the finite difference computation of Hessf

% OUTPUTS:
% Hessfx = n-by-n matrix corresponding to the approximation of the Hessian 
% of f in x.


Hessfx = zeros(length(x), length(x));

for j=1:length(x)
    % compute the elements (j, j) of the Heassian (diagonal)
    xh_plus = x;
    xh_minus = x;
    xh_plus(j) = xh_plus(j) + h;
    xh_minus(j) = xh_minus(j) - h;
    Hessfx(j,j) = (f(xh_plus) - 2*f(x) + f(xh_minus))/(h^2);
    for i=(j+1):length(x)  % value of cols such that we are considering only the U part
        % compute the element (i, j) of the Hessian
        xh_plus_ij = x;
        xh_plus_ij([i, j]) = xh_plus_ij([i, j]) + h;
        xh_plus_i = x;
        xh_plus_i(i) = xh_plus_i(i) + h;
        xh_plus_j = x;
        xh_plus_j(j) = xh_plus_j(j) + h;
        Hessfx(i,j) = (f(xh_plus_ij) - ...
            f(xh_plus_i) - f(xh_plus_j) + f(x))/(h^2);
                
        Hessfx(j,i)=Hessfx(i,j);
    end
end


end