function y=ODEBEP(a,f,y0,t);
% ODEBEP solve ODE on a given time mesh with parallelized Backward Euler
%   y=ODEBEP(a,f,y0,t); solves the ordinary differential equation
%   y'+a*y=f(t) with y(0)=y0 using Backward Euler on the time
%   mesh t=[0,t1,...,tN] by diagonalization, which could be performed in
%   parallel. f here is a function. 

y(1)=y0;
N=length(t)-1;
dt=t(2:end)-t(1:end-1);
B=spdiags([-1./[dt(2:end) 1]' 1./dt'],[-1 0],N,N);
[V,E]=eig(full(B));                    % compute decomposition numerically
fv=f(t(2:end))';
fv(1)=fv(1)+1/(t(2)-t(1))*y0;
y(2:N+1)=(V*((E+a*speye(N))\(V\fv)))';