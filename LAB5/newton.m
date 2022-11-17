function [xk,fk,gradfk_norm,k,xseq] = newton(x0,f,gradf,Hessf,kmax,tolgrad)
%NEWTON Summary of this function goes here
%   Detailed explanation goes here
    xseq = [x0];

    xk = x0;
    k = 0;
    gradfk_norm = norm(gradf(xk));

    while k < kmax && gradfk_norm > tolgrad
        Hessf_k = Hessf(xk);
        gradf_k = gradf(xk);
        %Hessf_k*pk = -gradf_k);
        pk = Hessf_k\(-gradf_k);

        x0 = xk;

        % Update values
        xk = x0 + pk;
        fk = f(xk);
        gradfk_norm = norm(gradf(xk));
        k = k + 1;
        
        xseq = [xseq, xk];
    end
end

