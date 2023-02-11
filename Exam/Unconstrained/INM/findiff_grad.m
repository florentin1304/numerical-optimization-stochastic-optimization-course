function [gradfx] = findiff_grad(f, x, h, type)

% INPUTS:
% f = function handle that describes a function;
% x = vector of size nx1;
% h = the h used for the finite difference computation of gradf,
% approximation step;
% type = 'fw' or 'c';


% OUTPUTS:
% gradfx = column vector (same size of x) corresponding to the approximation
% of the gradient of f in x

% initialize gradfx
gradfx = zeros(size(x));

switch type
    case 'fw'
        % for each i from 1 to n, f_xi ~ (f (x + h*ei) - f(x))/h
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;  % x+h*ei
            gradfx(i) = (f(xh) - f(x))/ h;
            
        end
    case 'c'
        for i=1:length(x)
            xh_plus = x;
            xh_minus = x;
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h;
            gradfx(i) = (f(xh_plus) - f(xh_minus))/(2 * h);
            
        end
    otherwise % repeat the 'fw' case
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x))/h;
        end
end



end

