function [gf] = gradf5(x)
    n = length(x);
    p = 7/3;

    df_pred = @ (xpp, xp, xi) p*abs((3-2*xp)*xp-xpp-xi+1)^(p-1)*abs((3-2*xp)*xp-xpp-xi+1)/((3-2*xp)*xp-xpp-xi+1)*(-1);
    df_curr = @ (xp, xi, xs) p*abs((3-2*xi)*xi-xp-xs+1)^(p-1)*abs((3-2*xi)*xi-xp-xs+1)/((3-2*xi)*xi-xp-xs+1)*(3-4*xi);
    df_succ = @ (xi, xs, xss) p*abs((3-2*xs)*xs-xi-xss+1)^(p-1)*abs((3-2*xs)*xs-xi-xss+1)/((3-2*xs)*xs-xi-xss+1)*(-1);
    gf = zeros(n, 1);

    for i=1:n
        x_i = x(i);

        % i == 1 and i == 2 (i starts from 1 and not 0, problems with x_pre and x_pp)
        if i == 1
            x_pre = 0;
            df_p = 0;
        elseif i == 2
            x_pre = x(i-1);
            x_pp = 0;
            df_p = df_pred(x_pp, x_pre, x_i);
        else
            x_pre = x(i-1);
            x_pp = x(i-2);
            df_p = df_pred(x_pp, x_pre, x_i);
        end

        % i == n and i == n-1  (problems with x_succ and x_ss)
        if i == n
            x_succ = 0;
            df_s = 0;
        elseif i == n-1
            x_succ = x(i+1);
            x_ss = 0;
            df_s = df_succ(x_i, x_succ, x_ss);
        else
            x_succ = x(i+1);
            x_ss = x(i+2);
            df_s = df_succ(x_i, x_succ, x_ss);
        end

        df_c = df_curr(x_pre, x_i, x_succ);
        gf(i) = df_p  + df_c + df_s;
    end

    
end
