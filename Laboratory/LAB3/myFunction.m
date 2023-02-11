function [y,dydx1] = myFunction(x1)

y = sin((pi/4)*norm(x1)^2);
[dydx1] = dlgradient(y,x1);

end