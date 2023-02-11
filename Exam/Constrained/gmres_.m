function [w, gradL_norm, eq_norm] = gmres_(Q, c, A, b, tol, maxit, n, K)
    d = [-c; b];
    Z = zeros(K, K);
    mat_K = [Q A'; A Z];
    w = gmres(mat_K, d, [], tol, maxit);
    gradL_norm = norm(Q * w(1:n,:) + c + A' * w(n+1:end, :));
    eq_norm = norm(A * w(1:n,:) - b);
end