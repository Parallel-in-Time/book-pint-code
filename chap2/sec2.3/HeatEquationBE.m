function u=HeatEquationBE(f,tspan,xspan,u0,gl,gr);
% HEATEQUATIONBE solve heat equation with Backward Euler
%   u=HeatEquationBE(f,Tspan,xspan,u0,gl,gr); solves the heat equation
%   with forcing function f on the domain xspan and initial condition
%   in the vector u0 whose length determines the number of spatial
%   gridpoints, and Dirichlet boundary conditions in the vectors gl
%   and gr whose length determines the number of grid points in time
%   using backward Euler on the time interval tspan. The corner points
%   are taken from the initial condition, not the boundary
%   condition. The output matrix u contains the solution in time
%   and space.

J=length(u0)-1; N=length(gl)-1;                   % determine grid
dt=(tspan(2)-tspan(1))/N;                         % mesh parameters
dx=(xspan(2)-xspan(1))/J;
x=(xspan(1):dx:xspan(2));
t=(tspan(1):dt:tspan(2))';
u(1,:)=u0;                                        % store initial and
u(2:N+1,1)=gl(2:end); u(2:N+1,J+1)=gr(2:end);     % boundary data
e=ones(J-1,1);
A=spdiags([e -2*e e],[-1 0 1],J-1,J-1)/dx^2;      % set up Laplacian
for n=1:N,
  rhs=u(n,2:end-1)+dt*feval(f,x(2:end-1),t(n+1)); % source function
  rhs(1)=rhs(1)+dt/dx^2*gl(n+1);                  % boundary conditions
  rhs(end)=rhs(end)+dt/dx^2*gr(n+1);
  u(n+1,2:J)=((speye(J-1)-dt*A)\rhs')';           % matrix solve
end;
