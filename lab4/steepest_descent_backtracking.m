function [xk,fk,gradfk_norm,k,xseq,btseq] = steepest_descent_backtracking(x0,f, ...
                                gradf,alpha0,kmax,tolgrad,c1,rho,btmax) %also new inputs
%STEEPEST_DESCENT_BACKTRACKING 
%Serve quando la funzione è molto rumorosa (tipo funzione banana)
%Oppure quando il minimo è 'disponibile' per un passo brevissimo

    xseq = zeros(length(x0), kmax);
    btseq = zeros(btmax);

    xk = x0;
    k = 0;
    gradfk_norm = norm(gradf(xk));

    while k < kmax && gradfk_norm < tolgrad
        pk = -gradf(xk);
        x0 = xk;
        
        bt = 0;
        alphak = alpha0;
        xnew = x0 + alphak*pk; 
        fnew = f(xnew);
        while bt < btmax && fnew > f(x0) + c1*alphak*(-pk')*pk
            %update alpha
            alphak = rho*alphak;
            xnew = x0 + pk*alphak;
            fnew = f(xnew);
            bt = bt + 1;
        end
        
        xk = xnew;
        fk = f(xk);
        gradfk_norm = norm(gradf(xk));
        k = k + 1;
        
        xseq(:, k) = xk;
        btseq(k) = bt;
    end



    xseq = xseq(1:k);
    btseq = btseq(1:k);
end

