function u=WaveEquation(f,c,tspan,xspan,u0,u0t,gl,gr)
% WAVEEQUATION solves the wave equation with centered finite differences
%   u=WaveEquation(f,c,Tspan,xspan,u0,u0t,gl,gr); solves the wave
%   equation with wave speed c and forcing function f on the domain
%   xspan and initial conditions in the vector u0 and u0t whose length
%   determines the number of spatial gridpoints, and Dirichlet
%   boundary conditions in the vectors gl and gr whose length
%   determines the number of grid points in time. The corner points
%   are taken from the boundary condition, not the initial condition,
%   to avoid inacuracies due to the implementation of the derivative
%   initial condition. The output matrix u contains the solution in
%   time and space.

J=length(u0)-1; N=length(gl)-1;                     % determine grid
dt=(tspan(2)-tspan(1))/N; dx=(xspan(2)-xspan(1))/J; % mesh parameters
x=(xspan(1):dx:xspan(2)); t=(tspan(1):dt:tspan(2))';
u(1,:)=u0;                                          % store initial data
u(2,:)=u0+dt*u0t;                                   % derivative by Taylor
u(1:N+1,1)=gl(1:end); u(1:N+1,J+1)=gr(1:end);       % store boundary data
e=ones(J-1,1);
for n=2:N,
  u(n+1,2:J)=2*u(n,2:J)-u(n-1,2:J)+c^2*dt^2/dx^2*(u(n,3:J+1)...
    -2*u(n,2:J)+u(n,1:J-1))+dt^2*feval(f,x(2:end-1),t(n));
end;