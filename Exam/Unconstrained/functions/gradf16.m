function [gf] = gradf16(x)
    n = length(x);
    gf = zeros(n, 1);

    for i=1:n
        x_i = x(i);
        if i == n
           gf(i) = -(i-1)*cos(x_i)+i*sin(x_i);
        else
           gf(i) = -(i-1)*cos(x_i)+i*sin(x_i)+(i+1)*cos(x_i);
        end
    end    
end
