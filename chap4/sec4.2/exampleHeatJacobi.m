l=7; J=2^l-1;                                  % mesh points in space
e=ones(J,1); dx=1/(J+1); x=(0:dx:1)';          % spatial mesh
L=1/dx^2*spdiags([e -2*e e],[-1 0 1],J,J);     % discrete Laplacian
T=5; N=J; dt=T/N; t=(0:dt:T);                  % time mesh
u0=@(x) zeros(size(x));                        % initial condition
gl=@(t) zeros(size(t));                        % boundary conditions
gr=@(t) zeros(size(t));
f=@(x,t) x.^4.*(1-x).^4+10*sin(8*t);           % source function
for n=1:N                                      % compute rhs b for reuse
  b(:,n)=dt*feval(f,x(2:end-1),t(n+1));        % source function
  b(1,n)=b(1,n)+dt/dx^2*gl(t(n+1));            % boundary conditions
  b(end,n)=b(end,n)+dt/dx^2*gr(t(n+1));
end;
u(:,1)=u0(x);                                  % set initial and boundary
u(1,1:N+1)=gl(t); u(end,1:N+1)=gr(t);          % conditions
A=speye(size(L))-dt*L;                         % time stepping matrix
for n=1:N                                      % compute exact solution
  u(2:end-1,n+1)=A\(u(2:end-1,n)+b(:,n));      % exact BE
end
uBE=u;                                         % keep exact BE solution
D=diag(diag(A));
NU=[5000 1000 100 5]
for l=1:length(NU)
  nu=NU(l)                                     % use nu Jacobi steps
  for n=1:N
    v=u(2:end-1,n);                            % use initial guess from
    for j=1:nu                                 % previous time step
      v=v+D\(u(2:end-1,n)+b(:,n)-A*v);         % for Jacobi iteration
    end;                                       % with no damping
    u(2:end-1,n+1)=v;
  end
  mesh(x,t,u'); xlabel('x'); ylabel('t');
  title(['Approximation after nu=' num2str(nu) ' Jacobi steps']);
  axis([0 1 0 5 -1 1]); pause;
end