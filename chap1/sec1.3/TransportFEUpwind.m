function [u,x,t]=TransportFEUpwind(f,u0,ua,a,b,J,T,N);
% TRANSPORTFEUPWIND solves the 1d transport equation 
%   [u,x,t]=TransportFEUpwind(f,u0,ua,a,b,J,T,N); solves the 1d
%   transport equation u_t+u_x=f in (a,b)x(0,T) with right hand side
%   function f and given intial condition u0, given boundary
%   condition ua at a using upwind finite differences in space with
%   J+1 mesh points and forward Euler in time up to T using N
%   meshpoints in time.

dx=(b-a)/J;
dt=T/N;
x=(a:dx:b)';
t=(0:dt:T);
u(1,1:N+1)=ua(t);
u(1:J+1,1)=u0(x);
for n=1:N,
  u(2:J+1,n+1)=u(2:J+1,n)-dt/dx*(u(2:J+1,n)-u(1:J,n))+dt*f(x(2:J+1),n*dt);
end;