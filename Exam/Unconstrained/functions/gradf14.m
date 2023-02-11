function [gradfx] = gradf14(x)
    %F14 Summary 
    %   Detailed explanation 

    n = length(x);
    h = 1/(n+1);
    gradfx = zeros(size(x));
    df_prec =@ (xpp, xp, xc, i) 2*(2*xp - xpp - xc + (1/2)*(h^2)*((xp+i*h + 1)^3) ) * (-1);
    df_curr =@ (xp, xc, xs, i) 2*(2*xc - xp - xs + (1/2)*(h^2)*((xc+i*h + 1)^3) ) * (2 + (3/2)*( (xc+ i*h +1)^2 ));
    df_succ =@ (xc, xs, xss, i)   2*(2*xs - xc - xss + (1/2)*(h^2)*((xs+i*h + 1)^3) ) * (-1);
    
    for i=1:1:n
        xc = x(i);

        % Control the edge cases for i==1 and i==2
        if i == 1
            xp = 0;
            df_p = 0;
        elseif i == 2
            xp = x(i-1);
            xpp = 0;
            df_p = df_prec(xpp, xp, xc, i);
        else
            xp = x(i-1);
            xpp = x(i-2);
            df_p = df_prec(xpp, xp, xc, i);
        end
        
        % Control the edge cases for i==n and i==n-1
        if i == n
            xs = 0;
            df_s = 0;
        elseif i == n-1
            xs = x(i+1);
            xss = 0;
            df_s = df_succ(xc,xs,xss,i);
        else
            xs = x(i+1);
            xss = x(i+2);
            df_s = df_succ(xc,xs,xss,i);
        end
        
        % Compute the current
        df_c = df_curr(xp,xc,xs,i);

        gradfx(i) = df_p + df_c + df_s;
    end

end


