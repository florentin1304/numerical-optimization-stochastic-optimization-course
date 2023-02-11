function [gradfx] = findiff_grad(f,x,h,type)
%FINDIFF_GRAD Approximation of f's gradient in x

    n = length(x);
    E = eye(n,n);
    gradfx = zeros(size(x));
    
    if strcmpi(type,"FW")
        fx = f(x);
        for i=1:n
            gradfx(i) = ( f(x+E(:,i)*h) - fx) / h;
        end

    elseif strcmpi(type, "C")
        for i=1:n
            gradfx(i) = ( f(x+E(:,i)*h) - f(x-E(:,i)*h)) / (2*h);
        end

    else

        error('findiff_grad type must be FW (forward) or C (centered)')
    end
end

