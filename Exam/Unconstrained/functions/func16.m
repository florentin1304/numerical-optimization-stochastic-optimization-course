function [y] = func16(x)
    
    n = length(x);
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

        y = y + i*((1-(cos(x_i)))+(sin(x_pred))-(sin(x_succ)));
    end

end