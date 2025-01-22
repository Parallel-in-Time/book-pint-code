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

Jc=(J+1)/2-1;                                  % coarse grid points
P=sparse(J,Jc);                                % prolongation by
for j=1:Jc                                     % interpolation
  P(2*j,j)=1; P(2*j-1,j)=0.5; P(2*j+1,j)=0.5;
end;
R=0.5*P';                                      % restriction by transpose
Lc=R*L*P;                                      % coarse matrix by Galerkin
Ac=speye(size(Lc))-dt*Lc;                      % coarsening in space only
xc=x(1:2:end);                                 % coarse spatial mesh
                                        % random initial guess
u(2:end-1,2:end)=rand(J,N);                    % with correct ic and bc
nu=5; al=0.5; K=10;
for k=1:K
  for n=1:N                                    % presmooting
    v=u(2:end-1,n+1);                          % use initial guess from
    for j=1:nu                                 % previous time step
      v=v+al*(D\(u(2:end-1,n)+b(:,n)-A*v));    % for Jacobi iteration
    end;                                       % with damping alpha
    u(2:end-1,n+1)=v;
  end
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after ' num2str(nu) ' presmoothing steps, iteration k=' num2str(k)]);
  pause
  for n=1:N                                    % compute residual
    r(:,n)=u(2:end-1,n)+b(:,n)-A*u(2:end-1,n+1);
  end;
  rc=R*r;                                      % restrict residual in space
  uc=zeros(Jc+2,1);                            % zero ic for correction
  for n=1:N                                    % coarse correction
    uc(2:end-1,n+1)=Ac\(uc(2:end-1,n)+rc(:,n));% by exact BE
  end
  u(2:end-1,:)=u(2:end-1,:)+P*uc(2:end-1,:);   % add coarse correction
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after coarse correction, iteration k=' num2str(k)]);
  pause
  for n=1:N                                    % postsmooting
    v=u(2:end-1,n+1);                          % use initial guess from
    for j=1:nu                                 % previous time step
      v=v+al*(D\(u(2:end-1,n)+b(:,n)-A*v));    % for Jacobi iteration
    end;                                       % with damping alpha
    u(2:end-1,n+1)=v;
  end
  err(k)=max(max(abs(uBE-u)));
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after ' num2str(nu) ' postsmoothing steps, iteration k=' num2str(k)]);
  pause
end

Nc=(N+1)/2-1;                                  % coarse grid points in time
Pc=sparse(N,Nc);                               % prolongation by interpolation
for j=1:Nc
  Pc(2*j,j)=1; Pc(2*j-1,j)=0.5; Pc(2*j+1,j)=0.5;
end;
Rc=0.5*Pc';                                    % restriction by transpose
Pc(end,end)=1;                                 % no zero bc in time at the end
Ac=speye(size(Lc))-2*dt*Lc;                    % coarsening in time also
u=uBE;                                         % random initial guess
u(2:end-1,2:end)=rand(J,N);                    % with correct ic and bc
al=0.5; nu=5; K=10;
for k=1:K
  for n=1:N                                    % presmooting
    v=u(2:end-1,n+1);                          % use initial guess from
    for j=1:nu                                 % previous time step
      v=v+al*(D\(u(2:end-1,n)+b(:,n)-A*v));    % for Jacobi iteration
    end;                                       % with damping alpha
    u(2:end-1,n+1)=v;
  end
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after ' num2str(nu) ' presmoothing steps, iteration k=' num2str(k)]);
  pause
  for n=1:N                                    % compute residual
    r(:,n)=u(2:end-1,n)+b(:,n)-A*u(2:end-1,n+1);
  end;
  rc=R*r;                                      % restrict residual in space
  rc=rc*Rc';                                   % restrict residual in time
  uc=zeros(Jc+2,1);                            % zero ic for correction
  for n=1:Nc                                   % coarse correction
    uc(2:end-1,n+1)=Ac\(uc(2:end-1,n)+rc(:,n));% exact BE using 2*dt
  end
  u(2:end-1,2:end)=u(2:end-1,2:end)+P*uc(2:end-1,2:end)*Pc';
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after coarse correction, iteration k=' num2str(k)]);
  pause
  for n=1:N                                    % postsmooting
    v=u(2:end-1,n+1);                          % use initial guess from
    for j=1:nu                                 % previous time step
      v=v+al*(D\(u(2:end-1,n)+b(:,n)-A*v));    % for Jacobi iteration
    end;                                       % with damping alpha
    u(2:end-1,n+1)=v;
  end
  errc(k)=max(max(abs(uBE-u)));
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after ' num2str(nu) ' postsmoothing steps, iteration k=' num2str(k)]);
  pause
end
