function [t,u]=ForwardEuler(f,tspan,u0,N);
% FORWARDEULER solves system of ODEs using the Forward Euler method
%   [t,u]=ForwardEuler(f,tspan,u0,N) solves du/dt=f(t,u) with initial
%   value u0 on the time interval tspan doing N steps of Forward
%   Euler. Returns the solution in time and space in the matrix u, and
%   also the corresponding time points in the column vector t.

dt=(tspan(2)-tspan(1))/N;
t=(tspan(1):dt:tspan(2))';   
u(1,:)=u0(:);                           % colon to make column vector
for n=1:N,
  u(n+1,:)=u(n,:)+dt*f(t(n),u(n,:));
end;
