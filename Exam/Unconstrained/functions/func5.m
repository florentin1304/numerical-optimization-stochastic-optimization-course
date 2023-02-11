function [y] = func5(x)
    %F5 Summary of this function goes here
    %   Detailed explanation goes here
    n = length(x);
    p = 7/3;
    y = 0;
    
    for i=1:1:n
        x_i = x(i);
        if i == 1
            x_pre = 0;
        else
            x_pre = x(i-1);
        end

        if i == n
            x_succ = 0;
        else
            x_succ = x(i+1);
        end

        y = y + abs( (3-2*x_i)*x_i - x_pre - x_succ +1 )^p;
    end

end