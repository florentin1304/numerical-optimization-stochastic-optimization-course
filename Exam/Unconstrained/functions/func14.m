function [y] = func14(x)
    %F14 Summary 
    %   Detailed explanation 

    n = length(x);
    h = 1/(n+1);
    y = 0;
    
    for i=1:1:n
        x_i = x(i);
        if i == 1
            x_pred = 0;
        else
            x_pred = x(i-1);
        end

        if i == n
            x_succ = 0;
        else
            x_succ = x(i+1);
        end

        y = y + ( 2*x_i - x_pred - x_succ + (h^2) * ((x_i + i*h + 1)^3)/2 )^2;
    end

end
