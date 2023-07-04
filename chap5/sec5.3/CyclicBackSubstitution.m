function x=CyclicBackSubstitution(A,f,xr);
% CYCLICBACKSUBSTITUTION backsubstitution after cyclic reduction
%   x=CyclicBackSubstitution(A,f,xr); performs the backsubstitution
%   in the bidiagional system Ax=b when the even unknowns are already
%   computed and given in xr.

n=length(f);                               % must be even 
id=(1:2:n-1)';                             % odd indices
x(id+1,1)=xr;                              % insert known solution values
d=diag(A); dm=diag(A,-1);                  % extract needed diagonal entries
x(id)=(f(id+1)-d(id+1).*xr)./dm(id);       % compute odd solution values 
