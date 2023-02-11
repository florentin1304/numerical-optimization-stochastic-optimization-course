function [y] = func13(x)
    %F13 Summary 
    %   Detailed explanation 
    n = length(x);
    k = n/2;
    y = 0;
    
    for j=1:1:k
        i = 2*j;
        x_i = x(i);
        x_pred = x(i-1);     

        y = y + (((x_pred)^2).^(x_i^2 + 1)+ (x_i^2)^(x_pred^2 + 1));
    end

end
