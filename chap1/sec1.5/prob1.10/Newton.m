function [eps, xk] = Newton(f, fA, x0, nIter)
% Define variables
eps = [];
xk = x0;

% Compute first right hand side and error
b = f(xk);
eps = [eps, norm(b)];

% Newton loop
for k = 1:nIter-1
    A = fA(xk);
    xkDiff = A\b;
    xk = xk - xkDiff;
    b = f(xk);
    eps = [eps, norm(b)];
end

end

