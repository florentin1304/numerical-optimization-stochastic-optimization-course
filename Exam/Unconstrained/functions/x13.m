function [x] = x13(n)
    x = ones(n,1);
    for i=1:1:n
        if mod(i, 2)==1
            x(i) = -1;
        elseif mod(i, 2)==0
            x(i) = 1;
        end
    end
end