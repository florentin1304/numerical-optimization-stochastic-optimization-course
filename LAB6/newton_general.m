function [xk,fk,gradfk_norm,k,xseq,btseq] = newton_general(x0,f,gradf,Hessf,kmax,tolgrad,c1,rho,btmax,verbose,FDgrad,FDhess,h)
%NEWTON Newton method with approximated gradient and hessian
    xseq = [x0];
    btseq = [];
        
    switch upper(FDgrad)
        case 'FW'
            gradf = @(x) findiff_grad(f, x, h, 'fw');
            
        case 'C'
            gradf = @(x) findiff_grad(f, x, h, 'c');
            
        otherwise
            
    end
    
    switch upper(FDhess)
        case 'FW'
            Hessf = @(x) findiff_Hess(f, x, sqrt(h));

        otherwise
            
    end


    xk = x0;
    k = 0;
    gradfk_norm = norm(gradf(xk));

    while k < kmax && gradfk_norm > tolgrad
        
        gradf_k = gradf(xk);
        Hessf_k = Hessf(xk);
     
        %Hessf_k*pk = -gradf_k;
        pk = Hessf_k\(-gradf_k);
        
        if verbose == true && not((gradf_k' * pk) < 0)
            warning('At k=%d -> pk is not descent vector!', k);
        end
        x0 = xk;
        
        bt = 0;
        alphak = 1;
        xnew = xk + alphak*pk;
        fnew = f(xnew);
        
        while bt < btmax && fnew > f(x0) + c1*alphak*(gradf(x0)')*pk
            alphak = rho*alphak;
            xnew = x0 + pk*alphak;
            fnew = f(xnew);
            bt = bt + 1;
        end

        xk = xnew;
        fk = fnew;
        gradfk_norm = norm(gradf(xk));
        k = k + 1;
        
        xseq = [xseq, xk];
        btseq = [btseq, bt];
    end
end

