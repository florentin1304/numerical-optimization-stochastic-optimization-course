%%%%%%%%%%%%%%%%%%%%% PARAMETERS
n = 1000;
alpha = 4;

%%%%%%%%%%%%%%%%%%%%% BUILD A SPARSE MATRIX
% choose the diagonals you want to use (0 = the square diagonal)
d = [-1, 0 ,1]; 

% give diagonal values
u = ones(n, 1); 
B = [-u, alpha*u, -u];

% create the diagonal with spdiags
A = spdiags(B,d,n,n);
x = ones(n, 1);
b = A*x;

%%%%%%%%%%%%%%%%%%%%% APPLY THE GRADIENT METHOD 
% stop method tollerance/max_iterations
tol = 10e-6;
kmax = 10e3;

% create initials 
x0 = zeros(n, 1); %(all zeros) 
r0 = b - A*x0; %(initial residuals)

% iterate
for i = 0:1:kmax
    zk = A*r0; % compute and save A*rk (faster)
    ak = ( (r0')*r0 )/( (r0')*zk ); % compute the step
    xk = x0 + ak*r0; % compute new 
    % rk = b - A*xk; 
    rk = r0 - ak*zk;
    
    % relative norm vs tollerance
    if norm(rk-r0, 2)/norm(rk, 2) < toll
       break 
    end
    
    % update values
    x0 = xk; % here not usefull
    r0 = rk;
end