function u=TransportBE(f,a,tspan,xspan,u0,N);
% TRANSPORTBE solve periodic transport equation with Backward Euler
%   u=TransportBE(f,a,Tspan,xspan,u0,N); solves the transport equation
%   u_t+au_x=f on the domain [0,L] and initial condition in the vector
%   u0 whose length determines the number of spatial gridpoints, and
%   periodic boundary conditions with N grid points in time using
%   backward Euler and upwind on the time interval tspan. The output
%   matrix u contains the solution in space and time.

J=length(u0)-1;                                  % determine grid
u(1,:)=u0;                                       % store initial data
dt=(tspan(2)-tspan(1))/N;                        % mesh parameters
dx=(xspan(2)-xspan(1))/J;
x=(xspan(1):dx:xspan(2));
t=(tspan(1):dt:tspan(2))';
e=ones(J,1);
A=spdiags([-e e],[-1 1],J,J)/(2*dx);              % set up periodic centered
A(1,end)=-1/(2*dx); A(end,1)=1/(2*dx);            % first derivative in x
for n=1:N,
  rhs=u(n,1:end-1)+dt*feval(f,x(1:end-1),t(n+1)); % add source term
  u(n+1,1:J)=((speye(J)+a*dt*A)\rhs')';           % matrix solve
  u(n+1,J+1)=u(n+1,1);                            % add periodic condition
end;