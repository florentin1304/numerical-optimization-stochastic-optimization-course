function [x] = x14(n)
    
    h = 1/(n+1);
    x = ones(n, 1);
    for i= 1:1:n
        x(i) = (i*h*(1-i*h));
    end
end


