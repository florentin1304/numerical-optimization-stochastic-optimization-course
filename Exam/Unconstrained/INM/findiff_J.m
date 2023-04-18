function [JFx] = findiff_J(F, x, h, type)
%
% function [JFx] = findiff_J(F, x, h, type)
%
% Function that approximate the Jacobian of F in x (column vector) with the
% finite difference (forward/centered) method.
%
% INPUTS:
% F = function handle that describes a function R^n->R^m;
% x = n-dimensional column vector;
% h = the h used for the finite difference computation of gradf
% type = 'fw' or 'c' for choosing the forward/centered finite difference
% computation of the gradient.
%
% OUTPUTS:
% JFx = matrix m-by-n corresponding to the approximation
% of the Jacobian of F in x.


JFx = zeros(length(F(x)), length(x));


switch type
    case 'fw'
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            JFx(:, i) = (F(xh) - F(x)) / h;
        end
    case 'c'
        for i=1:length(x)
            xh_plus = x;
            xh_minus = x;
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h;
            JFx(:, i) = (F(xh_plus) - F(xh_minus)) / (2 * h);
        end
    otherwise % same 'fw' case
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            JFx(:, i) = (F(xh) - F(x)) / h;
        end
end



end

