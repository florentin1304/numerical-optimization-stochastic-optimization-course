function [gradfx] = gradf13(x)

    n = length(x);
    gradfx = zeros(size(x));
    
%     df =@(xi,xoth) ( (2*(xoth^2) + 2) * ( xi^(2*(xoth^2)+1) )   +  (xoth^((2*(xi^2))+2)) * (4*xi) * log(xoth) );
    df =@(xi,xoth) 2*xi*(((xoth^2) + 1)*((xi^2)^(xoth^2))   +  ((xoth^2)^((xi^2)+1))*log(xoth^2)     );

    if mod(n,2) == 1
        n = n-1;
    end

    for i=1:1:n
        xi = x(i);

        if mod(i,2) == 0 
            xoth = x(i-1);
        else
            xoth = x(i+1);
        end

        gradfx(i) = df(xi,xoth);

    end
end

