function [B,g]=CyclicReduction(A,f);
% CYCLICREDUCTION performs a cyclic reduction for a bidiagonal system
%   [B,g]=CyclicReduction(A,f); performs a cyclic reduction for a lower
%   bidiagonal system Ax=b of even size.

n=length(f);                   % must be even 
id=1:2:n-1;

keyboard
d=diag(A);
dm=diag(A,-1);
dn=-dm(id+1)/d(id(2:end))*dm(id+2)
B=spdiags([d(id+1) dm],[-1 0],n/2,n/2);
