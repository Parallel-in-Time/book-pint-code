function [B, g] = CyclicReduction(A, f)
% CYCLICREDUCTION performs a cyclic reduction for a bidiagonal system
%   [B, g] = CyclicReduction(A, f); performs a cyclic reduction for a lower
%   bidiagonal system Ax = f of even size.

n = length(f);                   % must be even
iOdd = 1:2:n-1;
iEven = iOdd+1;

d = diag(A);                     % Main diagonal
dm = diag(A, -1);                % Sub-diagonal

dn = - dm(iOdd(2:end)) ./ d(iOdd(2:end)) .* dm(iEven(1:end-1));
B = spdiags([d(iEven), [dn; 0]], [0, -1], n/2, n/2);
g = f(iEven) - dm(iOdd) ./ d(iOdd) .* f(iOdd);
